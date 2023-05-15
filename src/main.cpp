#include "sentinel1_packet_decode.h"
#include "sentinel1.h"

int main() {
    Sentinel1PacketDecode sentinel1PacketDecode;
    sentinel1PacketDecode.ReadSARParam("C:/S1A_S3_RAW__0SDH_20220710T213600_20220710T213625_044043_0541DB_56CE/S1A_S3_RAW__0SDH_20220710T213600_20220710T213625_044043_0541DB_56CE.SAFE/s1a-s3-raw-s-hh-20220710t213600-20220710t213625-044043-0541db.dat");


    return 0;
}