echo "build index: fora foraplus foratp fpsn taupush";
./approx_dnppr -f amazon -alg fora -build 1 -sample 10 -random 1;
./approx_dnppr -f amazon -alg foraplus -build 1 -sample 10 -random 1;
./approx_dnppr -f amazon -alg foratp -build 1 -sample 10 -random 1;
./approx_dnppr -f amazon -alg fpsn -build 1 -sample 10 -random 1;
./approx_dnppr -f amazon -alg taupush -build 1 -sample 10 -random 1;

