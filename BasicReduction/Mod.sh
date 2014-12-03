ls -1 *.py > temp.txt

cat temp.txt | while read name
do
	sed -i "1s/^/#! \/usr\/bin\/python\n/" $name
done
