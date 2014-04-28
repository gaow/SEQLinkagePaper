import os, sys
import fnmatch
PROJ = "seqlink"
LOCAL = "tigerwang@tigerwang.org:./html/data/pages/software/{0}".format(PROJ)
PUB = "wanggao@bioinformatics.org:./data/{0}/pages".format(PROJ)
class Publisher:
    def __init__(self, baseurl):
        self.url = baseurl
        self.extsource = ".notes"
        self.extpage = ".txt"

    def compilePage(self, pages):
        pages = [os.path.splitext(x)[0] + self.extsource for x in pages]
        for page in pages:
            self.compile(page)

    def uploadPage(self, pages):
        pages = [os.path.splitext(x)[0] + self.extpage for x in pages]
        self.upload(self.url, " ".join(pages))

    def compilePages(self, path):
        for root, dirnames, filenames in os.walk(path):
            for filename in fnmatch.filter(filenames, "*" + self.extsource):
                self.compile(os.path.join(root, filename))

    def uploadPages(self, path):
        pages = []
        for root, dirnames, filenames in os.walk(path):
            for filename in fnmatch.filter(filenames, "*" + self.extpage):
                pages.append(os.path.join(root, filename))
        self.upload(self.url, " ".join(pages))

    def upload(self, url, pages):
        print("Uploading {0} ...".format(pages))
        os.system("scp {1} {0}".format(url, pages))

    def compile(self, page):
        print("Compiling {0} ...".format(page))
        os.system("tigernotes dokuwiki {0} --lite --showall {1}".format(page, '--toc' if page != 'start.notes' else ''))

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in ['local', 'pub']:
        sys.exit("USAGE: {0} {{local/pub}} [page [page ...]]".format(sys.argv[0]))
    if sys.argv[1] == "local":
        target = LOCAL
    else:
        target = PUB
    print("Publishing to {0}".format(target))
    p = Publisher(target)
    if len(sys.argv) == 2:
        p.compilePages(".")
        p.uploadPages(".")
    else:
        p.compilePage(sys.argv[2:])
        p.uploadPage(sys.argv[2:])
