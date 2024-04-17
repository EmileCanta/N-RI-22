for i in {1..13}
do
root -l -b -q "Sorter.C+(\"/Users/cantacuzene/data/n-ri-22/runs/raw_runs/TETRA23_RUN121.root\",\"/Users/cantacuzene/data/n-ri-22/runs/sorted_runs/RUN121_$i.root\",3000,6500,500,$i)"
done