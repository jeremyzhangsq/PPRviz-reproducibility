echo "twego plotting";
python2.7 benchmark.py --data 0 --repeat 1 --mode plot --algo pprvizs;
python2.7 benchmark.py --data 0 --repeat 1 --mode plot --algo fa2;
python2.7 benchmark.py --data 0 --repeat 1 --mode plot --algo ll;
python2.7 benchmark.py --data 0 --repeat 1 --mode plot --algo fr;
python2.7 benchmark.py --data 0 --repeat 1 --mode plot --algo mds;
python3 benchmark.py --data 0 --repeat 1 --mode plot --algo pivot_comp;
python2.7 benchmark.py --data 0 --repeat 1 --mode plot --algo pivot;
python3 benchmark.py --data 0 --repeat 1 --mode plot --algo simrank_comp;
python2.7 benchmark.py --data 0 --repeat 1 --mode plot --algo simrank;
python3 benchmark.py --data 0 --repeat 1 --algo maxent;
python2.7 benchmark.py --data 0 --repeat 1 --mode plot --algo gf;
python2.7 benchmark.py --data 0 --repeat 1 --mode plot --algo sdne;
python2.7 benchmark.py --data 0 --repeat 1 --mode plot --algo le;
python2.7 benchmark.py --data 0 --repeat 1 --mode plot --algo lle;
python2.7 benchmark.py --data 0 --repeat 1 --mode plot --algo node2vec;
echo "fbego plotting";
python2.7 benchmark.py --data 1 --repeat 1 --mode plot --algo pprvizs;
python2.7 benchmark.py --data 1 --repeat 1 --mode plot --algo fa2;
python2.7 benchmark.py --data 1 --repeat 1 --mode plot --algo ll;
python2.7 benchmark.py --data 1 --repeat 1 --mode plot --algo fr;
python2.7 benchmark.py --data 1 --repeat 1 --mode plot --algo mds;
python3 benchmark.py --data 1 --repeat 1 --mode plot --algo pivot_comp;
python2.7 benchmark.py --data 1 --repeat 1 --mode plot --algo pivot;
python3 benchmark.py --data 1 --repeat 1 --mode plot --algo simrank_comp;
python2.7 benchmark.py --data 1 --repeat 1 --mode plot --algo simrank;
python3 benchmark.py --data 1 --repeat 1 --algo maxent;
python2.7 benchmark.py --data 1 --repeat 1 --mode plot --algo gf;
python2.7 benchmark.py --data 1 --repeat 1 --mode plot --algo sdne;
python2.7 benchmark.py --data 1 --repeat 1 --mode plot --algo le;
python2.7 benchmark.py --data 1 --repeat 1 --mode plot --algo lle;
python2.7 benchmark.py --data 1 --repeat 1 --mode plot --algo node2vec;



