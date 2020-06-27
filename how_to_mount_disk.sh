### DD physique Data1 ###

cd /media/jbogoin/
sudo mkdir -p Data1
sudo mount /dev/sdb1 Data1

# Pour démonter
#sudo umount Data1

### Accès au serveur SortieSeqpilot ###

# Necessite le paquet cifs-utils
#sudo apt install cifs-utils

cd /media/jbogoin/
sudo mkdir -p SortieSeqPilot
sudo mount -t cifs -o username=miseq,password=illumina \
    //10.160.204.112/SeqData/SortieSeqPilot /media/jbogoin/SortieSeqPilot

### Accès au serveur fedseq ###
cd /media/jbogoin/
sudo mkdir -p fedseq
sudo mount -t cifs -o username=4006720,password=Abeille*99 //srv-psl/fedseq /media/jbogoin/fedseq

