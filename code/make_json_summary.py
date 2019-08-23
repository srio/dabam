if __name__ == '__main__':
    from dabam import make_json_summary
    make_json_summary()

    # # dump summary
    #
    # out_list = dabam_summary_dictionary()
    #
    # # print(out,type(out))
    #
    # out_dict = {}
    #
    # for i,ilist in enumerate(out_list):
    #     print("analyzing entry: ",i+1)
    #     out_dict["entry_%03d"%ilist["entry"]] = ilist
    #
    # print(out_dict)
    #
    # j = json.dumps(out_dict, ensure_ascii=True, indent="    ")
    #
    # print(j)
    # f = open("dabam-summary.json", 'w')
    # f.write(j)
    # f.close()
    # print("File dabam-summary.json written to disk")
    #
    # # dm = dabam()
    # # dm.load(12)

    # h = load_json_summary("dabam-summary.json")
    # for key in h:
    #     print(key)
    #
    # out = dabam_summary_dictionary_from_json_indexation(surface="elliptical(detrended)", slp_err_from=None, slp_err_to=None, length_from=None, length_to=None)
    #
    # for ilist in out:
    #     print(ilist)