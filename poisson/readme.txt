Bluegene:
Login: edu-cmc-stud16-615-09@bluegene.hpc.cs.msu.ru
Password: 1123581321

How to login:
cd /mnt/d/Dropbox/Projects/CMC/Supercomp/2
ssh -i ../mgu edu-cmc-stud16-615-09@bluegene.hpc.cs.msu.ru

How to upload changes:
cd /mnt/d/Dropbox/Projects/CMC/Supercomp/2
scp -i ../mgu poisson/main.c edu-cmc-stud16-615-09@bluegene1.hpc.cs.msu.ru:dev/src/poisson/
scp -i ../mgu poisson/Makefile edu-cmc-stud16-615-09@bluegene1.hpc.cs.msu.ru:dev/src/poisson/

How to build:
cd /dev/src/poisson/
make main

How to run:
make llsubmit
llq



Lomonosov:
Login: tereshinvs_1854@lomonosov.parallel.ru
Password: 1123581321

How to login:
cd /mnt/d/Dropbox/Projects/CMC/Supercomp/2
ssh -i ../mgu tereshinvs_1854@lomonosov.parallel.ru

How to upload changes:
scp -i ../mgu poisson/main.c tereshinvs_1854@lomonosov.parallel.ru:dev/src/poisson/
scp -i ../mgu poisson/Makefile tereshinvs_1854@lomonosov.parallel.ru:dev/src/poisson/

How to build:
modules av
modules add impi
make main
cp main ~/_scratch
cd ~/_scratch
sbatch -n16 -N2 -p test impi ./main

How to run:
cp main ~/_scratch
cd ~/_scratch
sbatch -n16 -N2 -p test impi ./main
