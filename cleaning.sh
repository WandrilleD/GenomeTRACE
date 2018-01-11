
for i in `ls source/*.cpp`;
    do 
    HERE=$(grep -c $i required_cpp);
    #echo $i $HERE

    if [ $HERE -lt 1 ]; 
        then rm $i
    fi
done;

for i in `ls include/*.h`;
    do 
    HERE=$(grep -c $i required_h);
    #echo $i $HERE

    if [ $HERE -lt 1 ]; 
        then rm $i
    fi
done;