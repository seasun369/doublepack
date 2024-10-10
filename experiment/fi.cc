#include <iostream>
#include <chrono>
#include <thread>

#include "dp/circuits.h"
#include "misc.h"

#define DELIM std::cout << "========================================\n"

#define DEBUG true
#define THREAD true

#define PRINT(x) if (DEBUG) std::cout << x << "\n";

inline std::size_t ValidateN(const std::size_t n) {
  assert(n > 3 && (((n-1) % 4) == 0));
  return n;			
}

inline std::size_t ValidateId(const std::size_t id, const std::size_t n) {
  if ( id >= n )
	  throw std::invalid_argument("ID cannot be larger than number of parties");
  return id;
}

int main(int argc, char** argv) {
  if (argc < 6) {
    std::cout << "usage: " << argv[0] << " [N] [id] [size] [depth]\n";
    return 0;
  }

  std::size_t n = ValidateN(std::stoul(argv[1]));
  std::size_t t = (n - 1) / 2;
  std::size_t l = std::stoul(argv[2]);
  std::size_t id = ValidateId(std::stoul(argv[3]), n);
  std::size_t size = std::stoul(argv[4]);
  std::size_t depth = std::stoul(argv[5]);
  std::size_t width = size/depth;

  DELIM;
  std::cout << "Running benchmark with N " << n << ", size " <<
    size << ", width " << width << " and depth " << depth << "\n";
  DELIM;

  auto config = scl::NetworkConfig::Localhost(id, n);
  // std::cout << "Config:"
  //           << "\n";
  // for (const auto& party : config.Parties()) {
  //   std::cout << " -- " << party.id << ", " << party.hostname << ", "
  //             << party.port << "\n";
  // }

  std::cout << "Connecting ..."
            << "\n";
  auto network = scl::Network::Create(config);

  std::cout << "Done!\n";

  std::size_t batch_size = (t + 2)/2;
  std::size_t n_parties = t + 2*(batch_size - 1) + 1;
  // std::size_t n_clients = n_parties;

  dp::CircuitConfig circuit_config;
  circuit_config.n_parties = n_parties;
  circuit_config.inp_gates = std::vector<std::size_t>(n_parties, 0);
  circuit_config.inp_gates[0] = 2;
  circuit_config.out_gates = std::vector<std::size_t>(n_parties, 0);
  circuit_config.out_gates[0] = 2;
  circuit_config.width = width;
  circuit_config.depth = depth;
  circuit_config.m_size = batch_size;
  circuit_config.l_size = l;
  circuit_config.batch_size = l * batch_size;
    
  auto circuit = dp::Circuit::FromConfig(circuit_config);

  circuit.SetNetwork(std::make_shared<scl::Network>(network), id);

  circuit.GenCorrelator();
  circuit.SetThreshold(t);

  DELIM;
  std::cout << "Running function-independent preprocessing\n";
  
  START_TIMER(fi_prep);
  PRINT("fi_prep SEND");
  if (THREAD) {
    std::thread t_FIPrepSend( &dp::Circuit::FIPrepSend, &circuit );
    PRINT("fi_prep RECV");
    circuit.FIPrepRecv();
    t_FIPrepSend.join();
  } else {
    circuit.FIPrepSend();
    PRINT("fi_prep RECV");
    circuit.FIPrepRecv();
  }

  PRINT("fi_prod");
  if (THREAD) {
    std::thread t_GenProdPartiesSendP1( &dp::Circuit::GenProdPartiesSendP1, &circuit ); 
    std::thread t_GenProdP1ReceivesAndSends( &dp::Circuit::GenProdP1ReceivesAndSends, &circuit ); 
    circuit.GenProdPartiesReceive();
    t_GenProdPartiesSendP1.join();
    t_GenProdP1ReceivesAndSends.join();
  } else {
    circuit.GenProdPartiesSendP1();
    circuit.GenProdP1ReceivesAndSends();
    circuit.GenProdPartiesReceive();
  }

  STOP_TIMER(fi_prep);

    std::cout << "\nclosing the network ...\n";
  // technically not necessary as channels are destroyed when their dtor is
  // called.
  network.Close();

}