if [ $# -lt 1 ]; then
    echo "usage: bash $0 version path"
    exit 0
fi
set -e
rm source/__init__.py
svn up
python setup.py install --prefix=$2
rm -rf dist
python development/release.py --prefix=$2 -p $2 --version $1
rsync -avP --dirs dist/* wanggao@bioinformatics.org:./public_html/seqlink/download/release/$1/
