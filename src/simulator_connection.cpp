#include "simulator_connection.h"

#include "simulator_session.h"

#include <boost/beast/core.hpp>
#include <boost/beast/websocket.hpp>
#include <exception>
#include <iostream>
#include <utility>

void ServeSimulatorConnection(boost::asio::ip::tcp::socket socket) {
  namespace beast = boost::beast;
  namespace websocket = beast::websocket;
  websocket::stream<boost::asio::ip::tcp::socket> connection(std::move(socket));
  connection.read_message_max(64 * 1024);
  connection.accept();
  connection.text(true);
  std::cout << "Simulator connected; starting a fresh run." << std::endl;

  SimulatorSession session;
  for (;;) {
    beast::flat_buffer buffer;
    beast::error_code error;
    connection.read(buffer, error);
    if (error == websocket::error::closed) {
      return;
    }
    if (error) {
      throw beast::system_error(error);
    }
    if (!connection.got_text()) {
      connection.close(websocket::close_code::unknown_data);
      return;
    }

    std::optional<std::string> response;
    try {
      response = session.HandleMessage(beast::buffers_to_string(buffer.data()));
    } catch (const std::exception &exception) {
      std::cerr << "Rejected simulator message: " << exception.what() << '\n';
      response = "42[\"manual\",{}]";
    }
    if (response) {
      connection.write(boost::asio::buffer(*response));
    }
  }
}
