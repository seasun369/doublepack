#include <catch2/catch.hpp>
#include "scl/net/channel.h"
#include "scl/net/network.h"
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_p.h>

TEST_CASE("Specialized Send and Recv for ZZ_pE and ZZ_p", "[channel]") {
    auto networks = scl::Network::CreateFullInMemory(2);
    auto channel0 = networks[0].Party(1);
    auto channel1 = networks[1].Party(0);

    SECTION("ZZ_pE Send and Recv") {
        NTL::ZZ_p::init(NTL::ZZ(17)); // Set a prime modulus
        NTL::ZZ_pX P;
        NTL::SetCoeff(P, 0, 1);
        NTL::SetCoeff(P, 1, 1);
        NTL::ZZ_pE::init(P); // Set the polynomial for ZZ_pE

        NTL::ZZ_pE sent_value;
        random(sent_value);

        channel0->Send(sent_value);

        NTL::ZZ_pE received_value;
        channel1->Recv(received_value);

        REQUIRE(sent_value == received_value);
    }

    SECTION("ZZ_p Send and Recv") {
        NTL::ZZ_p::init(NTL::ZZ(23)); // Set a different prime modulus

        NTL::ZZ_p sent_value;
        random(sent_value);

        channel0->Send(sent_value);

        NTL::ZZ_p received_value;
        channel1->Recv(received_value);

        REQUIRE(sent_value == received_value);
    }
}
