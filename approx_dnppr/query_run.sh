echo "query: taupush fpsn foratp fora+ resacc fora powiter";
./approx_dnppr -f amazon -alg taupush -build 0 -sample 10 -random 1;
./approx_dnppr -f amazon -alg fpsn -build 0 -sample 10 -random 1;
./approx_dnppr -f amazon -alg foratp -build 0 -sample 10 -random 1;
./approx_dnppr -f amazon -alg foraplus -build 0 -sample 2 -random 1;
./approx_dnppr -f amazon -alg foraresacc -build 0 -sample 2 -random 1;
./approx_dnppr -f amazon -alg fora -build 0 -sample 10 -random 1;
./approx_dnppr -f amazon -alg powiter -build 0 -sample 10 -random 1;

