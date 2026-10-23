#include "simulator_connection.h"
#include "simulator_session.h"

#include <boost/asio/io_context.hpp>
#include <boost/beast/core.hpp>
#include <boost/beast/websocket.hpp>
#include <catch2/catch_test_macros.hpp>
#include <exception>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {
using nlohmann::json;
std::string Telemetry(const std::string &row) {
  return "42" + json::array({"telemetry", {{"sensor_measurement", row}}}).dump();
}

json Payload(const std::optional<std::string> &reply) {
  if (!reply || reply->compare(0, 2, "42") != 0) {
    throw std::runtime_error("Missing Socket.IO response");
  }
  const auto event = json::parse(reply->substr(2));
  if (event[0] != "estimate_marker") {
    throw std::runtime_error("Missing estimate_marker event");
  }
  return event[1];
}

// A real TCP/WebSocket handshake and message exchange, with no Unity dependency.
std::vector<std::string> Exchange(const std::vector<std::string> &messages) {
  namespace net = boost::asio;
  namespace beast = boost::beast;
  using net::ip::tcp;
  net::io_context context;
  tcp::acceptor acceptor(context, {net::ip::make_address("127.0.0.1"), 0});
  const auto endpoint = acceptor.local_endpoint();
  std::exception_ptr server_error;
  std::thread server([&] {
    try {
      ServeSimulatorConnection(acceptor.accept());
    } catch (...) {
      server_error = std::current_exception();
    }
  });
  beast::websocket::stream<tcp::socket> client(context);
  std::vector<std::string> replies;
  try {
    client.next_layer().connect(endpoint);
    client.handshake("127.0.0.1:" + std::to_string(endpoint.port()),
                     "/socket.io/?EIO=3&transport=websocket");
    for (const auto &message : messages) {
      client.write(net::buffer(message));
      beast::flat_buffer buffer;
      client.read(buffer);
      replies.push_back(beast::buffers_to_string(buffer.data()));
    }
    client.close(beast::websocket::close_code::normal);
  } catch (...) {
    boost::system::error_code ignored;
    client.next_layer().close(ignored);
    server.join();
    throw;
  }
  server.join();
  if (server_error) {
    std::rethrow_exception(server_error);
  }
  return replies;
}
} // namespace

TEST_CASE("Simulator accepts full and short sensor rows", "[simulator]") {
  SimulatorSession full;
  SimulatorSession short_rows;
  const auto first = full.HandleMessage(Telemetry("L 2 0 1000000 2 0 0 0 0 0"));
  REQUIRE(first == short_rows.HandleMessage(Telemetry("L 2 0 1000000 2 0 0 0")));
  const auto second = full.HandleMessage(Telemetry("R 2 0 0 1100000 2 0 0 0 0 0"));
  REQUIRE(second == short_rows.HandleMessage(Telemetry("R 2 0 0 1100000 2 0 0 0")));
  const auto result = Payload(second);
  REQUIRE(result.size() == 6);
  REQUIRE(result["estimate_x"] == 2.0);
  REQUIRE(result["estimate_y"] == 0.0);
  for (const auto *name : {"rmse_x", "rmse_y", "rmse_vx", "rmse_vy"}) {
    REQUIRE(result[name] == 0.0);
  }
}

TEST_CASE("Simulator rejects malformed telemetry without changing the run", "[simulator]") {
  SimulatorSession session;
  SimulatorSession control;
  const auto first = Telemetry("L 2 0 1000000 2 0 0 0");
  REQUIRE(session.HandleMessage(first) == control.HandleMessage(first));
  REQUIRE_THROWS(session.HandleMessage("42[broken"));
  REQUIRE_THROWS(session.HandleMessage("42[\"telemetry\",{}]"));
  REQUIRE_THROWS(session.HandleMessage(Telemetry("R -2 0 0 1100000 2 0 0 0")));
  const auto next = Telemetry("L 2.1 0 1200000 2.1 0 0 0");
  REQUIRE(session.HandleMessage(next) == control.HandleMessage(next));
}

TEST_CASE("Simulator resets filter and RMSE together when a run restarts", "[simulator]") {
  SimulatorSession session;
  session.HandleMessage(Telemetry("L 2 0 1000000 3 0 0 0"));
  session.HandleMessage(Telemetry("L 2.1 0 1100000 3 0 0 0"));
  SECTION("Timestamp rewinds") {
    const auto restarted = Payload(session.HandleMessage(Telemetry("L 5 1 1000000 5 1 0 0")));
    REQUIRE(restarted["estimate_x"] == 5.0);
    REQUIRE(restarted["rmse_x"] == 0.0);
  }
  SECTION("Explicit reset event") {
    REQUIRE(session.HandleMessage("42[\"reset\",{}]") == "42[\"manual\",{}]");
    const auto restarted = Payload(session.HandleMessage(Telemetry("L 5 1 1200000 5 1 0 0")));
    REQUIRE(restarted["estimate_x"] == 5.0);
    REQUIRE(restarted["rmse_x"] == 0.0);
  }
}

TEST_CASE("Simulator handles control messages and null telemetry", "[simulator]") {
  SimulatorSession session;
  REQUIRE(session.HandleMessage("2probe") == "3probe");
  REQUIRE(session.HandleMessage("40") == "40");
  REQUIRE(session.HandleMessage("42null") == "42[\"manual\",{}]");
  REQUIRE(session.HandleMessage("42[\"telemetry\",null]") == "42[\"manual\",{}]");
  REQUIRE_FALSE(session.HandleMessage("42[\"unrelated\",{}]").has_value());
  REQUIRE_FALSE(session.HandleMessage("unrelated").has_value());
}

TEST_CASE("Simulator serves WebSocket telemetry and fresh connections", "[simulator][network]") {
  const auto replies = Exchange(
      {"2", Telemetry("L 2 0 1000000 3 0 0 0"), "42[broken", Telemetry("R 2 0 0 1100000 2 0 0 0")});
  REQUIRE(replies.size() == 4);
  REQUIRE(replies[0] == "3");
  REQUIRE(Payload(replies[1])["rmse_x"] == 1.0);
  REQUIRE(replies[2] == "42[\"manual\",{}]");
  REQUIRE(Payload(replies[3])["estimate_x"] == 2.0);
  const auto fresh = Exchange({Telemetry("L 5 0 1000000 5 0 0 0")});
  REQUIRE(Payload(fresh[0])["estimate_x"] == 5.0);
  REQUIRE(Payload(fresh[0])["rmse_x"] == 0.0);
}
