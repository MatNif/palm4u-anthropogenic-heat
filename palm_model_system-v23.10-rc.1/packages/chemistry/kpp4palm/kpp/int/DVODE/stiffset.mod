V24 stiffset
16 stiffoptions.f90 S582 0
05/10/2011  17:18:32
enduse
D 33 21 9 1 3 12 0 0 0 0 0
 0 12 3 3 12 12
D 36 21 9 1 3 13 0 0 0 0 0
 0 13 3 3 13 13
D 39 21 9 1 3 15 0 0 1 0 0
 0 14 3 3 15 15
D 42 21 9 1 3 15 0 0 1 0 0
 0 14 3 3 15 15
D 45 21 9 1 3 17 0 0 1 0 0
 0 16 3 3 17 17
D 48 21 9 2 18 21 0 0 1 0 0
 0 19 3 3 20 20
 0 16 20 3 17 17
D 51 21 9 1 3 23 0 0 1 0 0
 0 22 3 3 23 23
D 54 21 9 2 24 27 0 0 1 0 0
 0 25 3 3 26 26
 0 22 26 3 23 23
D 57 21 9 1 3 0 0 0 0 1 0
 0 0 3 3 0 28
D 60 21 6 1 3 0 0 0 0 1 0
 0 0 3 3 0 29
D 63 21 6 1 3 0 0 0 0 1 0
 0 0 3 3 0 30
D 66 21 9 1 3 0 0 0 0 1 0
 0 0 3 3 0 31
D 69 21 9 1 3 12 0 0 0 0 0
 0 12 3 3 12 12
D 72 21 9 1 3 12 0 0 0 0 0
 0 12 3 3 12 12
D 75 21 9 1 3 12 0 0 0 0 0
 0 12 3 3 12 12
D 78 21 9 1 3 12 0 0 0 0 0
 0 12 3 3 12 12
D 81 21 9 1 3 12 0 0 0 0 0
 0 12 3 3 12 12
S 582 24 0 0 0 8 1 0 4668 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 33 0 0 0 0 0 0 stiffset
S 583 6 4 0 0 6 584 582 4677 4 0 0 0 0 0 0 0 0 0 594 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 iwt
S 584 6 4 0 0 6 585 582 4681 4 0 0 4 0 0 0 0 0 0 594 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 n
S 585 6 4 0 0 6 586 582 4683 4 0 0 8 0 0 0 0 0 0 594 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 id
S 586 6 4 0 0 6 587 582 4686 4 0 0 12 0 0 0 0 0 0 594 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 iid
S 587 6 4 0 0 6 588 582 4690 4 0 0 16 0 0 0 0 0 0 594 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 nfcn
S 588 6 4 0 0 6 589 582 4695 4 0 0 20 0 0 0 0 0 0 594 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 njac
S 589 6 4 0 0 6 1 582 4700 4 0 0 24 0 0 0 0 0 0 594 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 nlud
S 590 7 4 0 4 33 591 582 4705 800004 100 0 0 0 0 0 0 0 0 595 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 w
S 591 7 4 0 4 36 1 582 4707 800004 100 0 160 0 0 0 0 0 0 595 0 0 0 0 0 0 0 0 0 0 582 0 0 0 0 dy
S 592 3 0 0 0 6 1 1 0 0 0 0 0 0 0 0 20 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 593 3 0 0 0 6 0 1 0 0 0 0 0 0 0 0 400 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 594 11 0 0 0 8 1 582 4710 40800000 801000 0 28 0 0 583 589 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 stiffset$0
S 595 11 0 0 1 8 594 582 4721 40800000 801000 0 3360 0 0 590 591 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 stiffset$2
S 596 23 5 0 0 0 601 582 4732 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 derivs
S 597 6 3 0 0 6 1 596 4739 800004 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 neq
S 598 1 3 0 0 9 1 596 4743 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 t
S 599 7 3 0 0 39 1 596 4745 800204 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 y
S 600 7 3 0 0 42 1 596 4747 800204 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ydot
S 601 14 5 0 0 0 1 596 4732 200 400000 0 0 0 2 4 0 0 0 0 0 0 0 0 0 0 0 0 48 0 582 0 0 0 0 derivs
F 601 4 597 598 599 600
S 602 6 1 0 0 6 1 596 4752 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_14
S 603 23 5 0 0 0 611 582 4759 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 jacd
S 604 6 3 0 0 6 1 603 4739 800004 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 neq
S 605 1 3 0 0 9 1 603 4743 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 t
S 606 7 3 0 0 45 1 603 4745 800204 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 y
S 607 1 3 0 0 6 1 603 4764 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ml
S 608 1 3 0 0 6 1 603 4767 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 mu
S 609 7 3 0 0 48 1 603 4770 800204 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 pd
S 610 6 3 0 0 6 1 603 4773 800004 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nrowpd
S 611 14 5 0 0 0 1 603 4759 200 400000 0 0 0 7 7 0 0 0 0 0 0 0 0 0 0 0 0 58 0 582 0 0 0 0 jacd
F 611 7 604 605 606 607 608 609 610
S 612 6 1 0 0 6 1 603 4780 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_16
S 613 6 1 0 0 6 1 603 4787 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_18
S 614 6 1 0 0 6 1 603 4794 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_20
S 615 6 1 0 0 6 1 603 4801 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_22
S 616 23 5 0 0 0 624 582 4808 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 jacb
S 617 6 3 0 0 6 1 616 4739 800004 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 neq
S 618 1 3 0 0 9 1 616 4743 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 t
S 619 7 3 0 0 51 1 616 4745 800204 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 y
S 620 1 3 0 0 6 1 616 4764 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ml
S 621 1 3 0 0 6 1 616 4767 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 mu
S 622 7 3 0 0 54 1 616 4770 800204 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 pd
S 623 6 3 0 0 6 1 616 4773 800004 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nrowpd
S 624 14 5 0 0 0 1 616 4808 200 400000 0 0 0 15 7 0 0 0 0 0 0 0 0 0 0 0 0 73 0 582 0 0 0 0 jacb
F 624 7 617 618 619 620 621 622 623
S 625 6 1 0 0 6 1 616 4801 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_22
S 626 6 1 0 0 6 1 616 4813 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_24
S 627 6 1 0 0 6 1 616 4820 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_26
S 628 6 1 0 0 6 1 616 4827 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_28
S 629 23 5 0 0 0 637 582 4834 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 jacs
S 630 1 3 0 0 6 1 629 4739 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 neq
S 631 1 3 0 0 9 1 629 4743 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 t
S 632 7 3 0 0 57 1 629 4745 800104 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 y
S 633 7 3 0 0 60 1 629 4839 800104 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ia
S 634 7 3 0 0 63 1 629 4842 800104 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ja
S 635 1 3 0 0 6 1 629 4845 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nz
S 636 7 3 0 0 66 1 629 4848 800104 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 p
S 637 14 5 0 0 0 1 629 4834 100 400000 0 0 0 23 7 0 0 0 0 0 0 0 0 0 0 0 0 88 0 582 0 0 0 0 jacs
F 637 7 630 631 632 633 634 635 636
S 638 6 1 0 0 6 1 629 4850 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_0
S 639 6 1 0 0 6 1 629 4856 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_1
S 640 6 1 0 0 6 1 629 4862 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_2
S 641 6 1 0 0 6 1 629 4868 40800006 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_3
S 642 23 5 0 0 0 648 582 4874 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ivalu
S 643 1 3 0 0 9 1 642 4880 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 xstart
S 644 1 3 0 0 9 1 642 4887 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 xend
S 645 1 3 0 0 9 1 642 4892 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 hbegin
S 646 1 3 0 0 9 1 642 4899 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 hmax
S 647 7 3 0 0 69 1 642 4745 800004 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 y
S 648 14 5 0 0 0 1 642 4874 0 400000 0 0 0 31 5 0 0 0 0 0 0 0 0 0 0 0 0 118 0 582 0 0 0 0 ivalu
F 648 5 643 644 645 646 647
S 649 23 5 0 0 0 651 582 4904 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 evalu
S 650 7 3 0 0 72 1 649 4745 800004 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 y
S 651 14 5 0 0 0 1 649 4904 0 400000 0 0 0 37 1 0 0 0 0 0 0 0 0 0 0 0 0 609 0 582 0 0 0 0 evalu
F 651 1 650
S 652 23 5 0 0 0 656 582 4910 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 fcn
S 653 1 3 0 0 9 1 652 4914 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x
S 654 7 3 0 0 75 1 652 4745 800004 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 y
S 655 7 3 0 0 78 1 652 4916 800004 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 yp
S 656 14 5 0 0 0 1 652 4910 0 400000 0 0 0 39 3 0 0 0 0 0 0 0 0 0 0 0 0 889 0 582 0 0 0 0 fcn
F 656 3 653 654 655
S 657 23 5 0 0 0 660 582 4919 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 pderv
S 658 1 3 0 0 9 1 657 4914 4 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 x
S 659 7 3 0 0 81 1 657 4745 800004 3000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 y
S 660 14 5 0 0 0 1 657 4919 0 400000 0 0 0 43 2 0 0 0 0 0 0 0 0 0 0 0 0 1135 0 582 0 0 0 0 pderv
F 660 2 658 659
A 12 2 0 0 0 6 592 0 0 0 12 0 0 0 0 0 0 0 0 0
A 13 2 0 0 0 6 593 0 0 0 13 0 0 0 0 0 0 0 0 0
A 14 1 0 0 0 6 597 0 0 0 0 0 0 0 0 0 0 0 0 0
A 15 1 0 0 0 6 602 0 0 0 0 0 0 0 0 0 0 0 0 0
A 16 1 0 0 0 6 604 0 0 0 0 0 0 0 0 0 0 0 0 0
A 17 1 0 0 0 6 612 0 0 0 0 0 0 0 0 0 0 0 0 0
A 18 1 0 0 0 6 615 0 0 0 0 0 0 0 0 0 0 0 0 0
A 19 1 0 0 0 6 610 0 0 0 0 0 0 0 0 0 0 0 0 0
A 20 1 0 0 0 6 613 0 0 0 0 0 0 0 0 0 0 0 0 0
A 21 1 0 0 0 6 614 0 0 0 0 0 0 0 0 0 0 0 0 0
A 22 1 0 0 0 6 617 0 0 0 0 0 0 0 0 0 0 0 0 0
A 23 1 0 0 0 6 625 0 0 0 0 0 0 0 0 0 0 0 0 0
A 24 1 0 0 0 6 628 0 0 0 0 0 0 0 0 0 0 0 0 0
A 25 1 0 0 0 6 623 0 0 0 0 0 0 0 0 0 0 0 0 0
A 26 1 0 0 0 6 626 0 0 0 0 0 0 0 0 0 0 0 0 0
A 27 1 0 0 0 6 627 0 0 0 0 0 0 0 0 0 0 0 0 0
A 28 1 0 0 0 6 638 0 0 0 0 0 0 0 0 0 0 0 0 0
A 29 1 0 0 0 6 639 0 0 0 0 0 0 0 0 0 0 0 0 0
A 30 1 0 0 0 6 640 0 0 0 0 0 0 0 0 0 0 0 0 0
A 31 1 0 0 0 6 641 0 0 0 0 0 0 0 0 0 0 0 0 0
Z
Z
