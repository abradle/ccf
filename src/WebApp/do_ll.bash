while [[ $# > 1 ]]
do
key="$1"

case $key in
    -m|--mol2_prot)
    MOLTWO="$2"
    shift # past argument
    ;;
    -p|--products)
    PRODUCTS="$2"
    shift # past argument
    ;;
    -l|--ll_prot)
    LLPROT="$2"
    shift # past argument
    ;;
    -s|--smiles)
    SMILES="$2"
    shift # past argument
    ;;
    -c|--context)
    CONTEXT="$2"
    shift # past argument
    ;;
    -n|--name)
    NAME="$2"
    shift # past argument
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done
sudo docker pull abradle/lloommppaa
sudo docker run -it -v $PWD/data:/data --name $NAME abradle/lloommppaa /bin/bash -c 'source ~/.bashrc && service postgresql restart && cd /CHOC/src/WebApp && python djangorun.py syncdb --noinput && python manage.py --targ LL_DO --mol2_prot '$MOLTWO' --products '$PRODUCTS' --ll_prot '$LLPROT' --lloommppaa True --smiles "'$SMILES'" --context "'$CONTEXT'" > out.log'