pip install pysam --user
pip install biopython --user
pip install intervaltree --user
pip install jellyfish --user
pip install numpy --user
pip install networkx --user
tar -xf suffix_tree-2.1.tar.gz
cd suffix_tree-2.1
python setup.py install --user
cd ..
rm -rf suffix_tree-2.1
