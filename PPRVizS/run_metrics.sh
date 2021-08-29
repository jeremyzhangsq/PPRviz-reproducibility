echo "twego plotting";
python2.7 benchmark.py --data 0 --repeat 10 --mode metrics --algo pprvizs;
python2.7 benchmark.py --data 0 --repeat 10 --mode metrics --algo fa2;
python2.7 benchmark.py --data 0 --repeat 10 --mode metrics --algo ll;
python2.7 benchmark.py --data 0 --repeat 10 --mode metrics --algo fr;
python2.7 benchmark.py --data 0 --repeat 10 --mode metrics --algo mds;
python3 benchmark.py --data 0 --repeat 10 --mode metrics --algo pivot_comp;
python2.7 benchmark.py --data 0 --repeat 10 --mode metrics --algo pivot;
python3 benchmark.py --data 0 --repeat 10 --mode metrics --algo simrank_comp;
python2.7 benchmark.py --data 0 --repeat 10 --mode metrics --algo simrank;
python3 benchmark.py --data 0 --repeat 10 --algo maxent;
python2.7 benchmark.py --data 0 --repeat 10 --mode metrics --algo gf;
python2.7 benchmark.py --data 0 --repeat 10 --mode metrics --algo sdne;
python2.7 benchmark.py --data 0 --repeat 10 --mode metrics --algo le;
python2.7 benchmark.py --data 0 --repeat 10 --mode metrics --algo lle;
python2.7 benchmark.py --data 0 --repeat 10 --mode metrics --algo node2vec;
echo "fbego plotting";
python2.7 benchmark.py --data 1 --repeat 10 --mode metrics --algo pprvizs;
python2.7 benchmark.py --data 1 --repeat 10 --mode metrics --algo fa2;
python2.7 benchmark.py --data 1 --repeat 10 --mode metrics --algo ll;
python2.7 benchmark.py --data 1 --repeat 10 --mode metrics --algo fr;
python2.7 benchmark.py --data 1 --repeat 10 --mode metrics --algo mds;
python3 benchmark.py --data 1 --repeat 10 --mode metrics --algo pivot_comp;
python2.7 benchmark.py --data 1 --repeat 10 --mode metrics --algo pivot;
python3 benchmark.py --data 1 --repeat 10 --mode metrics --algo simrank_comp;
python2.7 benchmark.py --data 1 --repeat 10 --mode metrics --algo simrank;
python3 benchmark.py --data 1 --repeat 10 --algo maxent;
python2.7 benchmark.py --data 1 --repeat 10 --mode metrics --algo gf;
python2.7 benchmark.py --data 1 --repeat 10 --mode metrics --algo sdne;
python2.7 benchmark.py --data 1 --repeat 10 --mode metrics --algo le;
python2.7 benchmark.py --data 1 --repeat 10 --mode metrics --algo lle;
python2.7 benchmark.py --data 1 --repeat 10 --mode metrics --algo node2vec;



