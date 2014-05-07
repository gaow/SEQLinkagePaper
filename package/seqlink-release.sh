if [ $# -lt 1 ]; then
    echo "usage: bash $0 version path"
    exit 0
fi
set -e
rm source/__init__.py
svn up
rm -rf dist
python development/release.py --prefix=$2 -p $2 --version $1
ssh wanggao@bioinformatics.org "cp -a public_html/seqlink/download/.template public_html/seqlink/download/$1"
rsync -avP --dirs dist/* wanggao@bioinformatics.org:./public_html/seqlink/download/$1
