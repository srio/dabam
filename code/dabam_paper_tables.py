import dabam

if __name__ == '__main__':

    print(dabam.dabam_summary(nmax=25,latex=0))

    txt = dabam.dabam_summary(nmax=25,latex=1)
    output_file="mytable1.tex"
    f1 = open(output_file,'w')
    f1.write(txt)
    f1.close()
    print("File written to disk %s"%output_file)

    txt = dabam.dabam_summary(nmax=25,latex=2)
    output_file="mytable2.tex"
    f1 = open(output_file,'w')
    f1.write(txt)
    f1.close()
    print("File written to disk %s"%output_file)
