#ifndef SIMULATOR_CONNECTION_H_
#define SIMULATOR_CONNECTION_H_

#include <boost/asio/ip/tcp.hpp>

// Accepts the WebSocket handshake and serves one connection until disconnect.
// Every connection gets a fresh EKF and RMSE history.
void ServeSimulatorConnection(boost::asio::ip::tcp::socket socket);

#endif
