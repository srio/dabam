import matplotlib.pyplot as plt

fig, ax = plt.subplots(2, 2)


#
# surface type
#

if True:
    tt=[
        "       plane  ",
        "       plane  ",
        "   spherical  ",
        "  elliptical  ",
        "   spherical  ",
        "  elliptical  ",
        "       plane  ",
        "       plane  ",
        "       plane  ",
        "       plane  ",
        "       plane  ",
        "       plane  ",
        "   spherical  ",
        "   spherical  ",
        "    toroidal  ",
        "    toroidal  ",
        "    toroidal  ",
        " cylindrical  ",
        "  elliptical  ",
        "  elliptical  ",
        "  elliptical  ",
        " cylindrical  ",
        "       plane  ",
        "       plane  ",
        "       plane  ",
        "   spherical  ",
        "    Toroidal  ",
        "       Plane  ",
        "   Spherical  ",
        "       Plane  ",
        "       Plane  ",
        "       Plane  ",
        "       Plane  ",
        "       Plane  ",
        "       Plane  ",
        "       Plane  ",
        "   Spherical  ",
        "   Spherical  ",
        "elliptical",
        "elliptical",
        "   spherical  ",
        "elliptical",
        "toroid",
        "toroid",
        "    toroidal  ",
        "    toroidal  ",
        "    toroidal  ",
        "    toroidal  ",
        "    toroidal  ",
        "    toroidal  ",
        "       plane  ",
        "       plane  ",
        "       plane  ",
        "       plane  ",
        "toroid",
        "toroid",
        "toroid",
        "       Plane  ",
        "       Plane  ",
        "       Plane  ",
        "       Plane  ",
        "  Elliptical  ",
        "  Elliptical  ",
        "  Elliptical  ",
        "  Elliptical  ",
        "  Elliptical  ",
        "  Elliptical  ",
        "  Elliptical  ",
        "  Elliptical  ",
        "  Elliptical  ",
        "       Plane  ",
        "       Plane  ",
        "       Plane  ",
        "       Plane  ",
        "       Plane  ",
        "       Plane  ",
        "       Plane  ",
        "       Plane  ",
        " cylindrical  ",
        " cylindrical  ",
        "elliptical",
        "elliptical",
        "           1  ",
        "           1  ",
        "           1  ",
        "           1  ",
        "           1  ",
        "           1  ",
        "           1  ",
        "           1  ",
        "           1  ",
        "           1  ",
        "           1  ",
        "   spherical  ",
        "elliptical",
        "elliptical"]

    print(len(tt))

    count_type_plane = 0
    count_type_ellipse = 0
    count_type_toroid = 0
    count_type_sphere = 0
    count_type_cylinder = 0
    count_type_other = 0

    for tt1 in tt:
        print(tt1)
        if "lane" in tt1:
            count_type_plane +=1
            print("+plane")
        elif  "llip" in tt1:
            count_type_ellipse +=1
            print("+ellipse")
        elif  "oroid" in tt1:
            count_type_toroid +=1
            print("+toroid")
        elif  "pher" in tt1:
            count_type_sphere +=1
            print("+sphe")
        elif  "ylind" in tt1:
            count_type_cylinder +=1
            print("+cylinder")
        else:
            count_type_other += 1
            print("+other")

    print(count_type_plane, count_type_ellipse, count_type_toroid, count_type_sphere, count_type_cylinder, count_type_other)
    print(count_type_plane+count_type_ellipse+count_type_toroid+count_type_sphere+count_type_cylinder+count_type_other)

    labels = (
        "plane",
        "ellipse",
        "toroid",
        "sphere",
        "cylinder",
        "other")
    sizes = [count_type_plane, count_type_ellipse, count_type_toroid, count_type_sphere, count_type_cylinder, count_type_other]

    # fig, ax = plt.subplots()
    print(ax)
    ax[0,0].pie(sizes, labels=labels, autopct='%1.f%%')
    ax[0,0].set_title("Mirror curvature")
    # plt.show()


#
# length
#
if True:
    length = [
        1200.00,
        360.00 ,
        118.00 ,
        32.00  ,
        429.00 ,
        200.00 ,
        99.00  ,
        97.00  ,
        97.00  ,
        442.00 ,
        445.00 ,
        442.00 ,
        114.00 ,
        240.00 ,
        800.00 ,
        239.00 ,
        495.00 ,
        330.00 ,
        350.00 ,
        121.00 ,
        430.00 ,
        1130.00,
        127.00 ,
        774.00 ,
        200.00 ,
        948.00 ,
        289.00 ,
        439.00 ,
        185.00 ,
        480.00 ,
        299.00 ,
        299.00 ,
        390.00 ,
        239.00 ,
        239.00 ,
        239.00 ,
        1000.00,
        1000.00,
        1000.00,
        1000.00,
        700.00 ,
        900.00 ,
        598.00 ,
        598.00 ,
        1300.00,
        1300.00,
        1300.00,
        1200.00,
        1200.00,
        1200.00,
        190.00 ,
        190.00 ,
        140.00 ,
        140.00 ,
        900.00 ,
        600.00 ,
        700.00 ,
        220.00 ,
        220.00 ,
        220.00 ,
        220.00 ,
        90.00  ,
        47.00  ,
        240.00 ,
        240.00 ,
        480.00 ,
        480.00 ,
        480.00 ,
        240.00 ,
        240.00 ,
        968.00 ,
        968.00 ,
        968.00 ,
        968.00 ,
        972.00 ,
        972.00 ,
        972.00 ,
        972.00 ,
        269.00 ,
        270.00 ,
        170.00 ,
        200.00 ,
        500.00 ,
        294.00 ,
        444.00 ,
        294.00 ,
        444.00 ,
        294.00 ,
        444.00 ,
        234.00 ,
        849.00 ,
         80.00 ,
         80.00 ,
        120.00 ,
        59.00  ,
        25.00  ,
    ]


    print(len(length))

    count_length_150 = 0
    count_length_300 = 0
    count_length_600 = 0
    count_length_1200 = 0
    count_length_more = 0


    for tt1 in length:
        print(tt1)
        if tt1 <= 150:
            count_length_150 +=1
            print("+150")
        elif tt1 <= 300:
            count_length_300 +=1
            print("+300")
        elif tt1 <= 600:
            count_length_600 +=1
            print("+600")
        elif tt1 <= 1200:
            count_length_1200 +=1
            print("+1200")
        else:
            count_length_more +=1
            print("+more")



    print(count_length_150 ,
        count_length_300 ,
        count_length_600 ,
        count_length_1200,
        count_length_more)
    print(count_length_150 +
        count_length_300 +
        count_length_600 +
        count_length_1200+
        count_length_more)

    import matplotlib.pyplot as plt
    labels = (
        "0-150 mm",
        "150-300 mm",
        "300-600 mm",
        "600-1200 mm",
        "more")
    sizes = [count_length_150, count_length_300, count_length_600, count_length_1200, count_length_more]

    # fig, ax = plt.subplots()
    ax[0,1].pie(sizes, labels=labels, autopct='%1.f%%')
    ax[0,1].set_title("Mirror length")
    # plt.show()


#
# slope error
#
if True:
    slope = [
        0.49,
        0.15,
        0.17,
        0.13,
        0.84,
        1.01,
        0.58,
        0.57,
        0.44,
        0.23,
        0.23,
        0.32,
        1.30,
        1.15,
        2.13,
        1.85,
        2.27,
        0.12,
        0.14,
        0.45,
        0.43,
        5.32,
        0.18,
        0.20,
        0.07,
        0.86,
        1.86,
        0.95,
        0.60,
        0.89,
        0.56,
        1.33,
        0.29,
        0.44,
        0.45,
        0.48,
        0.32,
        0.53,
        0.34,
        0.35,
        0.44,
        1.78,
        0.92,
        1.95,
        0.78,
        0.75,
        0.86,
        0.80,
        1.50,
        0.49,
        0.38,
        0.61,
        0.34,
        0.38,
        1.60,
        0.53,
        0.74,
        0.17,
        0.19,
        0.23,
        0.19,
        0.51,
        0.91,
        0.75,
        0.69,
        1.77,
        1.53,
        1.52,
        0.45,
        0.58,
        0.26,
        0.30,
        0.20,
        0.23,
        0.27,
        0.23,
        0.27,
        0.19,
        0.16,
        0.16,
        0.59,
        0.31,
        0.15,
        0.18,
        0.33,
        0.18,
        0.16,
        0.17,
        0.13,
        0.07,
        0.62,
        0.99,
        0.25,
        0.16,
        0.36,
        0.40,
    ]

    print(len(slope))

    count_slope_p25  = 0
    count_slope_p5   = 0
    count_slope_1    = 0
    count_slope_2    = 0
    count_slope_more = 0


    for tt1 in slope:
        print(tt1)
        if tt1 <= 0.25:
            count_slope_p25 +=1
            print("+p25")
        elif tt1 <= 0.5:
            count_slope_p5 +=1
            print("+p5")
        elif tt1 <= 1:
            count_slope_1 +=1
            print("+1")
        elif tt1 <= 2:
            count_slope_2 +=1
            print("+1200")
        else:
            count_slope_more +=1
            print("+more")



    print(
    count_slope_p25 ,
    count_slope_p5  ,
    count_slope_1   ,
    count_slope_2   ,
    count_slope_more,
    )
    print(
    count_slope_p25 +
    count_slope_p5  +
    count_slope_1   +
    count_slope_2   +
    count_slope_more
          )

    import matplotlib.pyplot as plt
    labels = (
        r"0-0.25 $\mu$rad",
        r"0.25-0.5 $\mu$rad",
        r"0.5-1 $\mu$rad",
        r"1-2 $\mu$rad",
        r"more")
    sizes = [count_slope_p25, count_slope_p5, count_slope_1, count_slope_2, count_slope_more]

    # fig, ax = plt.subplots()
    ax[1,0].pie(sizes, labels=labels, autopct='%1.f%%')
    ax[1,0].set_title("slopes profile RMS")
    # plt.show()


#
# height rms
#
if True:
    heights = [
        43.85 ,
        4.46 ,
        1.56 ,
        0.22 ,
        31.90 ,
        2.98 ,
        0.51 ,
        4.05   ,
        3.52   ,
        6.20   ,
        6.59   ,
        2.86   ,
        8.97   ,
        23.28  ,
        170.87 ,
        54.78  ,
        61.54  ,
        1.83 ,
        7.25 ,
        3.40 ,
        6.05 ,
        400.82 ,
        2.17 ,
        5.35 ,
        1.47   ,
        94.52 ,
        41.98 ,
        51.73 ,
        5.22 ,
        37.00 ,
        14.93 ,
        61.28 ,
        6.80 ,
        11.66 ,
        6.93 ,
        7.74 ,
        11.95 ,
        44.99 ,
        20.41 ,
        22.64 ,
        16.34  ,
        174.91 ,
        16.99 ,
        25.98 ,
        101.37 ,
        89.86  ,
        101.70 ,
        80.45  ,
        186.95 ,
        34.44  ,
        8.14   ,
        10.09  ,
        5.41   ,
        4.25   ,
        153.30 ,
        18.15 ,
        41.71 ,
        1.12 ,
        1.38 ,
        2.14 ,
        1.37 ,
        0.65 ,
        2.24 ,
        19.65 ,
        17.80 ,
        51.57 ,
        47.34 ,
        35.91 ,
        11.56  ,
        15.10  ,
        5.45 ,
        4.43 ,
        3.56 ,
        5.36 ,
        3.11 ,
        5.70 ,
        5.60 ,
        4.57 ,
        1.25   ,
        1.34   ,
        10.38 ,
        5.60 ,
        1.58 ,
        2.30 ,
        5.29 ,
        1.65 ,
        6.18 ,
        2.47 ,
        5.16 ,
        1.07 ,
        7.04 ,
        8.32 ,
        2.17 ,
        2.07   ,
        0.48   ,
        0.35   ,
    ]
    print(len(heights))

    count_heights_1    = 0
    count_heights_3    = 0
    count_heights_6    = 0
    count_heights_12   = 0
    count_heights_more = 0


    for tt1 in heights:
        print(tt1)
        if tt1 <= 1:
            count_heights_1 +=1
            print("+1")
        elif tt1 <= 3:
            count_heights_3 +=1
            print("+3")
        elif tt1 <= 6:
            count_heights_6 +=1
            print("+6")
        elif tt1 <= 12:
            count_heights_12 +=1
            print("+12")
        else:
            count_heights_more +=1
            print("+more")



    print(
    count_heights_1   ,
    count_heights_3   ,
    count_heights_6   ,
    count_heights_12  ,
    count_heights_more,
    )
    print(
    count_heights_1   +
    count_heights_3   +
    count_heights_6   +
    count_heights_12  +
    count_heights_more
          )

    import matplotlib.pyplot as plt
    labels = (
        "0-1 nm",
        "1-3 nm",
        "3-6 nm",
        "6-12 nm",
        "more")
    sizes = [count_heights_1, count_heights_3, count_heights_6, count_heights_12, count_heights_more]

    # fig, ax = plt.subplots()
    ax[1,1].pie(sizes, labels=labels, autopct='%1.f%%')
    ax[1,1].set_title("height profile RMS")

    # plt.show()

plt.show()