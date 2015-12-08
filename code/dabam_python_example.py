import dabam

if __name__ == '__main__':


    #
    # get summary table
    #

    text = dabam.dabam_summary()
    print(text)


    out = dabam.dabam_summary_dictionary()
    for i in range(len(out)):
        print("i=:%d, entry:%d"%(i,(out[i])["entry"]))

    print(out[13])

    #
    # load a given entry (=14)
    #
    dm = dabam.dabam()
    dm.set_input_silent(True)
    dm.load(14)


    info = dm.info_profiles()
    print(info)

    #
    # make a bunh of plots
    #
    dm.set_input_plot("heights slopes psd_h csd_h acf_h histo_h")
    dm.plot()

    # you can do plots by accessing data, ex: plt.plot(1e3*dm.y,1e6*dm.zHeights)
    #     dm.metadata              # metadata
    #     dm.rawdata               # raw datafile
    #     dm.y                     # abscissa along the mirror
    #     dm.zSlopesUndetrended    # undetrended slope profile
    #     dm.zSlopes               # detrended slope profile
    #     dm.zHeightsUndetrended   # undetrended heights profile
    #     dm.zHeights              # detrended heights profile
    #     dm.coeffs                # information on detrending (polynomial coeffs)
    #     dm.f                     # frequency of Power Spectral Density
    #     dm.psdHeights            # Power Spectral Density of Heights profile
    #     dm.psdSlopes             # Power Spectral Density of slopes profile
    #     dm.csdHeights            # Antiderivative of PDF of Heights profile
    #     dm.csdSlopes             # Antiderivative of PDF of Slopes profile
    #     dm.histoSlopes           # to store slopes histogram
    #     dm.histoHeights          # to store heights histogram
    #     dm.momentsSlopes         # to store moments of the slopes profile
    #     dm.momentsHeights        # to store moments of the heights profile
    #     dm.powerlaw              {"hgt_pendent":None, "hgt_shift":None, "slp_pendent":None, "slp_shift":None,


