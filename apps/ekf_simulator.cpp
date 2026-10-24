#include "simulator_connection.h"

#include <boost/asio/io_context.hpp>
#include <charconv>
#include <exception>
#include <iostream>
#include <string_view>
#include <system_error>
#include <utility>

int main(int argc, char *argv[]) {
  unsigned short port = 4567;
  if (argc > 2) {
    std::cerr << "Usage: ekf_simulator [port]\n";
    return 1;
  }
  if (argc == 2) {
    const std::string_view argument(argv[1]);
    const auto result = std::from_chars(argument.data(), argument.data() + argument.size(), port);
    if (result.ec != std::errc{} || result.ptr != argument.data() + argument.size() || port == 0) {
      std::cerr << "Port must be an integer from 1 to 65535\n";
      return 1;
    }
  }
  try {
    using boost::asio::ip::tcp;
    boost::asio::io_context context;
    tcp::acceptor acceptor(context);
    acceptor.open(tcp::v4());
    acceptor.set_option(tcp::acceptor::reuse_address(true));
    acceptor.bind({boost::asio::ip::make_address("127.0.0.1"), port});
    acceptor.listen();
    std::cout << "Listening on ws://127.0.0.1:" << port
              << " — start the Ubuntu Term 2 simulator. Ctrl+C to stop." << std::endl;
    for (;;) {
      auto socket = acceptor.accept();
      try {
        ServeSimulatorConnection(std::move(socket));
      } catch (const std::exception &error) {
        std::cerr << "Simulator connection ended: " << error.what() << '\n';
      }
      std::cout << "Waiting for simulator connection..." << std::endl;
    }
  } catch (const std::exception &error) {
    std::cerr << "Simulator server failed: " << error.what() << '\n';
    return 1;
  }
}
