import os, sys
from subprocess import Popen, PIPE

def plotdemo(title, genotype = None):
  rdat = '''
    suppressMessages(library('kinship2'))
    data(sample.ped)'''
  docs = '''
    ## fam1 <- sample.ped[sample.ped$ped==1,]
    ## ped1 <- pedigree(fam1$id, fam1$father, fam1$mother,
    ## fam1$sex, fam1$affected, fam1$avail)
    ## ## to see plot:
    ## plot.pedigree(ped1, align=FALSE)
    ## id1 <- pedigree.unrelated(ped1, avail=fam1$avail)
    ## id1
    ## ## some possible vectors
    ## ##[1] "110" "113" "133" "109"
    ## ##[1] "113" "118" "141" "109"
    ## ##[1] "113" "118" "140" "109"
    ## ##[1] "110" "113" "116" "109"
    ## ##[1] "113" "133" "141" "109"
    ## fam2 <- sample.ped[sample.ped$ped==2,]
    ## ped2 <- pedigree(fam2$id, fam2$father, fam2$mother,
    ## fam2$sex, fam2$affected, fam2$avail)
    ## ## to see plot:
    ## ## plot.pedigree(ped2, align=FALSE)
    ## id2 <- pedigree.unrelated(ped2, avail=fam2$avail)
    ##'''
  fam1 = '''
    dat <- sample.ped[sample.ped$ped==2,]
    '''
  nuc1 = '''
    dat <- sample.ped[sample.ped$ped==2,]
    dat <- dat[dat$id %in% c(201,202,205,206,207),]
    dat$affected[2] <- 0
    dat$avail <- 0
    '''
  prt = '''
    dat$father[which(dat$father == 0)] <- NA
    dat$mother[which(dat$mother == 0)] <- NA
    dat <- pedigree(as.character(dat$id), as.character(dat$father), as.character(dat$mother),
    dat$sex, dat$affected, dat$avail)
    png('{0}.png')
    plot.pedigree(dat, symbolsize = 1, cex = 1.5, align = T)
    dev.off()
    '''.format(title)
  if genotype:
    edit = '''
    genotype <- %s
    if (length(genotype) > length(dat$id)) {
    genotype <- genotype[c(1,length(dat$id))]
    } else {
    # genotype <- c(genotype, dat$id[c(length(genotype)+1, length(dat$id))])
    genotype <- c(genotype, rep(NULL, length(dat$id) - length(genotype)))
    }
    id <- dat$id
    #
    for (i in 1:length(dat$id))
    dat$id[i] <- ifelse(is.null(genotype[i]), dat$id[i], paste(dat$id[i], genotype[i], sep='%s'))
    #
    for (i in 1:length(dat$father)) {
    for (j in 1:length(id)) {
    if (dat$father[i] == id[j]) {
    dat$father[i] = dat$id[j]
    break
    }
    }
    }
    #
    for (i in 1:length(dat$mother)) {
    for (j in 1:length(id)) {
    if (dat$mother[i] == id[j]) {
    dat$mother[i] = dat$id[j]
    break
    }
    }
    }
    ''' % ('c({})'.format(', '.join(['"{}"'.format('\n'.join(i)) for i in genotype])), '\n')
  else:
    edit = ''
  return rdat + nuc1 + edit + prt

def run(genotype, title):
  Popen(["R", "--slave"], stdin=PIPE).communicate(plotdemo(title, genotype))

def montage(title, tile = None):
  if tile:
    os.system("montage *.*.png -geometry 960x960+1+1 -tile %s %s" % (tile, title))
  else:
    os.system("montage *.*.png -geometry 960x960+1+1 %s" % (title))
  os.system("rm *.*.png")
  
if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit()
    if sys.argv[1] == '1':
        run([('0?','0?'),('0?','0?'),('??','00'),('00','??'),('??','??')], '1.1')
        run([('AB','00','01'),('CD','01','00'),('??','00'),('00','??'),('??','??')], '1.2')
        run([('AB','00','01'), ('CD','01','00'), ('AD','01','00'), ('CB','00','01'), ('BD','01','10')], '1.3')
        run([('12',), ('13',), ('13',), ('12',), ('23',)], '1.4')
        run([('AB','00','01'), ('CD','01','00'), ('DA','10','00'), ('CB','00','01'), ('BD','01','10')], '1.5')
        run([('12',), ('13',), ('31',), ('12',), ('23',)], '1.6')
        montage('Figure1.png', '2x3')
    if sys.argv[1] == '2':
        # run([('00','0?'),('01','0?'),('??','??'),('??','??'),('??','??')], '2.1')
        run([('AB','00','01','00','??'),('CD','01','00','01','??'),
             ('AD','01','00','??','??'),('CB','00','01','??','??'),('BD','01','10','??','??')], '2.1')

        montage('Figure2.png')
