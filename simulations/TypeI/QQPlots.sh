source $1
grep $gene $null | cut -f2 > /tmp/1.txt
grep $gene $recomb | cut -f2 > /tmp/2.txt
grep $gene $ldmissing | cut -f2 > /tmp/3.txt
echo """ 
library('ggplot2')
gghist <- function(dat, stat_foo = NULL, vlines = NULL, xname='x', verbose = F) {
        average <- round(mean(dat[,xname], na.rm=T),4)
        stdev <- round(sd(dat[,xname], na.rm=T),4)
        #med <- round(median(dat[,xname], na.rm=T), 4)
        r1 <- round(min(dat[,xname], na.rm=T),4)
        r2 <- round(max(dat[,xname], na.rm=T),4)
        # convert dat obj from numeric to data.frame
        myplot <- ggplot(dat) +
                        aes_string(x = xname) +
                        geom_histogram(aes_string(y = '..density..', fill = '..density..'), color = 'white', binwidth = (r2-r1)/30) +
                        scale_fill_gradient('bin mass', low = 'darkolivegreen3', high = colors()[552]) +
                        geom_rug(color = 'grey80') +
                        scale_x_continuous(name = paste('\\n', xname)) +
                        scale_y_continuous(name = 'Density\\n') +
                        theme_bw()
        if (!is.null(vlines)) myplot <- myplot + geom_vline(xintercept = vlines, color = '#9ACD32', linetype = 2)
        if (!is.null(stat_foo)) {
        myplot <- myplot + geom_line(stat = 'density', size = 0.5, linetype = 2, color = 'grey50') + stat_function(fun = stat_foo, color = 'red')
        plottitle <- 'Histogram & fitted vs. normal distribution density curve for '
        } else {
        plottitle <- 'Histogram for '
        }
        if (verbose) {
        myplot <- myplot + labs(title = paste(plottitle, xname, '\\n', 'mean = ', toString(average), '; ', 'stdev = ', toString(stdev), '; ', 'range = ', '[', toString(r1), ',', toString(r2), ']', '\\n', sep=''))
        } else {
        myplot <- myplot + labs(title = paste(plottitle, xname, '\\n', sep=''))
        }
        return(myplot)
}
ggdot <- function(dat, psize=2.5, color=NULL, xyz = NULL, title = '') {
    if (!is.null(xyz)) {
    colnames(dat) = xyz
    } else {
    xyz = colnames(dat)
    }
    myplot = ggplot(dat, aes_string(xyz[1], xyz[2])) +
        geom_point(size = psize, binwidth = range(dat[,1])/30) +
        xlab(paste('\\n', xyz[1], sep='')) + ylab(paste(xyz[2], '\\n', sep='')) +
        geom_abline(intercept = 0, slope = 1) +
	    theme(axis.title.x=element_text(size=10)) +
	    theme(axis.title.y=element_text(size=10)) +
	    theme_bw() +
        labs(title=paste(title,'\\n'))
    return(myplot) 
}
HLOD = scan('/tmp/1.txt')
len1 = length(HLOD)
d2 = scan('/tmp/2.txt')
len2 = length(d2)
d3 = scan('/tmp/3.txt')
len3 = length(d3)
d20 = sort(HLOD[1:min(len1, len2)])
d21 = sort(d2[1:min(len1, len2)])
d30 = sort(HLOD[1:min(len1, len3)])
d31 = sort(d3[1:min(len1, len3)])
print(length(which(HLOD>3.6)))
pdf('$gene.pdf')
print(gghist(data.frame(HLOD), xname='HLOD', vlines = c(3.6)))
print(ggdot(data.frame(d20,d21), xyz = c('Expected','Observed'), title = 'HLOD null distribution in the presence of recombination events'))
print(ggdot(data.frame(d30,d31), xyz = c('Expected','Observed'), title = 'HLOD null distribution in the presence of LD and missing data'))
dev.off()
""" | R --slave  
if [ "$(uname)" == "Darwin" ]; then
    prog='open'
else
    prog='google-chrome'
fi
$prog $gene.pdf 
