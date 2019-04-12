#g++ analyzer_stackRatio_order.cpp -I`root-config --incdir` `root-config --libs` -o plot
g++ analyzer_order_trig_stack_pas.cpp -I`root-config --incdir` `root-config --libs` -std=c++0x -o nopas

boolZ=0
if [ $2 = Z ]; then 
                boolZ=1
fi


./nopas 0 $1 $boolZ  #analysis

