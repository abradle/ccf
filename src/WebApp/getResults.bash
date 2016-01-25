sudo docker run -it -v $PWD/data:/data --name $NAME abradle/lloommppaa /bin/bash -c 'source ~/.bashrc && service postgresql restart && cd /CHOC/src/WebApp && python djangorun.py syncdb --noinput && python manage.py --targ LL_DO --mol2_prot '$MOLTWO' --products '$PRODUCTS' --ll_prot '$LLPROT' --lloommppaa True --smiles "'$SMILES'" --context "'$CONTEXT'" > out.log'
abradley@cocoa:~$ cat getResults.bash 
sudo docker commit $1 abradle/lloommppaa:$1
mkdir $1
cp new.py $1
echo $PWD/$1
sudo docker run -v $PWD/$1:/data -it abradle/lloommppaa:$1 /bin/bash -c "ls /data && service postgresql start && cp /data/new.py /CHOC/src/WebApp && source ~/.bashrc && cd /CHOC/src/WebApp && python new.py"