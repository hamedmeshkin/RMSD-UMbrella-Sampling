
# ################################################################ 
# ##### Reading From Files as initial parameter and so on ######## 
# ################################################################
set cmCoor0 {}
set fd [open input/Rindex.dat r]
set OK 1
while {1} {
    if {[gets $fd line] == -1} {
        break
    }   
    lappend cmCoor0 $line
    foreach aid $line {
        if {![info exists flag($aid)]} {
            addatom $aid
            set flag($aid) 1
        }
    }
}

set IndexList {}
set fd [open input/torsion_index.dat r]
set OK 1
while {1} {
    if {[gets $fd line] == -1} {
        break
    }
    if {[llength $line] != 4} {
        error "Every line must have 4 elements!"
    }
    lappend IndexList $line
    foreach aid $line {
        if {![info exists flag($aid)]} {
            addatom $aid
            set flag($aid) 1
        }
    }
}

set HelixList {}
set fd [open input/helix_idex.dat r]
set OK 1
while {1} {
    if {[gets $fd line] == -1} {
        break
    }
    if {[llength $line] != 12} {
        error "Every line must have 12 elements!"
    }
    lassign $line ind1 ind2 ind3 ind4 ind5 ind6 ind7 ind8 wr_HLX1 wr_HLX2 kw_HLX rang_HLX
    set lin "$ind1 $ind2 $ind3 $ind4 $ind5 $ind6 $ind7 $ind8"
    lappend HelixList $lin
    foreach aid $lin {
        if {![info exists flag($aid)]} {
            addatom $aid
            set flag($aid) 1
        }
    }
}
set kw_HLXTwic [expr $kw_HLX * 2.0]

set fd [open input/ref_win0.dat r]
set refList [read $fd]
##--------------------------------------------------------------------------##

# ############################################################################# 
# ##### Testin the Metropolis Criterion to Calculate Energy Difference ######## 
# #############################################################################
set nRef [llength $refList]
proc calcDiffEng {n} {  
    global refID phiD IndF angleFixer deg2Rad2 atmcrd nRef wr lable i1 i2 i3 i4 KwRM_half     
    global kw sec DelPhi ts qFile rang rangRM Energy distnce condR condT KwRM RmRef RMSD
    puts $qFile "$ts $RmRef $wr $sec"
    # If n is even, window 2*i exchanges with 2*i + 1;
    # If n is odd, window 2*i exchanges with 2*i - 1.
    if {($refID + $n) % 2 == 0} {
        set newRef [expr $refID + 1]
        if {$newRef >= $nRef} {
            return nan 
        }
    } else {
        set newRef [expr $refID - 1]
        if {$newRef < 0} {
            return nan
        }
    }

    IndF $newRef
    set DelPhi [angleFixer [expr ($phiD-$wr)]]   
    set DeltaRmsd [expr $RMSD - $RmRef]
    
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%% Bias Potential on the Center of Mass %%%%%%% 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if {$condR} {
        if {$lable=="MM"} {
            set EngA [expr  $KwRM_half*$DeltaRmsd*$DeltaRmsd]   
        } else {
            if {abs($DeltaRmsd) <= $rangRM} {
                set EngA  0.0
            } elseif  {$DeltaRmsd > 0.0} {
                set DeltaRmsd [expr ($DeltaRmsd-$rangRM)]
                set EngA [expr $KwRM_half*$DeltaRmsd*$DeltaRmsd] 
            } else {
                set DeltaRmsd [expr ($DeltaRmsd+$rangRM)]
                set EngA [expr $KwRM_half*$DeltaRmsd*$DeltaRmsd] 
            }
        }
        # Boundary
        if {abs($DelPhi) <= $rang} {
            set EngA_b   0.0 
        } elseif  {$DelPhi > 0.0} {
            set DelPhi  [expr ($DelPhi-$rang)]
            set EngA_b  [expr $kw*$DelPhi*$DelPhi*$deg2Rad2] 
        } else {
            set DelPhi [expr  ($DelPhi+$rang)]
            set EngA_b  [expr $kw*$DelPhi*$DelPhi*$deg2Rad2]
        }
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%% Bias Potential on the Dihedral Angle %%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    } elseif {$condT} {
        if {$lable=="TT"} {
            set EngA   [expr $kw*$DelPhi*$DelPhi*$deg2Rad2]
        } else {
            if {abs($DelPhi) <= $rang} {
                set EngA  0.0
            } elseif  {$DelPhi > 0.0} {
                set DelPhi [expr   ($DelPhi-$rang)]
                set EngA    [expr $kw*($DelPhi*$DelPhi)*$deg2Rad2]
            } else {
                set DelPhi [expr   ($DelPhi+$rang)]
                set EngA    [expr $kw*($DelPhi*$DelPhi)*$deg2Rad2]
            }
        }
        # Boundary     
        if {abs($DeltaRmsd) <= $rangRM} {
            set EngA_b   0.0 
        } elseif  {$DeltaRmsd > 0.0} {
            set DeltaRmsd   [expr ($DeltaRmsd-$rangRM)]
            set EngA_b [expr $KwRM_half * $DeltaRmsd * $DeltaRmsd] 
        } else {
            set DeltaRmsd   [expr ($DeltaRmsd+$rangRM)]
            set EngA_b [expr $KwRM_half * $DeltaRmsd * $DeltaRmsd] 
        }
    }
    return [expr $EngA+$EngA_b-$Energy]
} 
##-------------------------------------------------------------------##

# ##################################################################### 
# ##### Fix the Delta phi not to get a value bigger that |180| ######## 
# ##################################################################### 
proc angleFixer {DelPhi} {
    if {$DelPhi > 180} {
        set DelPhi [expr $DelPhi-360]
    } elseif {$DelPhi < -180} {
        set DelPhi [expr $DelPhi+360]
    }
    return $DelPhi
}
##----------------------------------------------------------------------------##

# ##############################################################################
# ##### Index of the Bias Potential and Assigning the related parameter ######## 
# ##############################################################################
proc IndF {refID} {
    global rang refList RmRef lable IndexList kw sec KwRM RMS
    global aid i1 i2 i3 i4 rangRM KwRM_half wr condR condT kw_twic
    set sec      $refID;  
    set lable    [lindex $refList $sec 0]
    set subsec   [lindex $refList $sec 1]
    set wr       [lindex $refList $sec 2]
    set kw       [lindex $refList $sec 3] ;# Kcal/mol/rad2
    set rang     [lindex $refList $sec 4]
    set RmRef    [lindex $refList $sec 5]
    set KwRM     [lindex $refList $sec 6] ;# Kcal/mol/A2
    set rangRM   [lindex $refList $sec 7]
    set RMS      [lindex $refList $sec 8]
    
    set aid      [lindex $IndexList  $subsec] 
    lassign $aid i1 i2 i3 i4   
    
    if {$lable=="TL"} {
        set wr [expr $wr-$rang]
    } elseif {$lable=="TR"} {
        set wr [expr $wr+$rang]
    }
    
    if {$lable=="ML"} {
        set RmRef [expr $RmRef-$rangRM]
    } elseif {$lable=="MR"} {
        set RmRef [expr $RmRef+$rangRM]
    }
    
    set condR [expr {$lable=="MM"} || {$lable=="ML"} || {$lable=="MR"}]
    set condT [expr {$lable=="TT"} || {$lable=="TL"} || {$lable=="TR"}]
    
    set kw_twic   [expr $kw * 2.0]
    set KwRM_half [expr $KwRM / 2.0]
}
##------------------------------------------------------------##

# ############################################### 
# ##### Angle & Aihedral Force On Torsion ####### 
# ###############################################
proc addForce {force kw0 T1 T2 T3 T4} {
    global  atmcrd ts
    set energy [expr $force*$force/$kw0]
    foreach {g1 g2 g3 g4} [dihedralgrad $atmcrd($T1) $atmcrd($T2) $atmcrd($T3) $atmcrd($T4)] {}
    addforce $T1 [vecscale $g1 $force]
    addforce $T2 [vecscale $g2 $force]
    addforce $T3 [vecscale $g3 $force]
    addforce $T4 [vecscale $g4 $force]
    return $energy    
} 
 
set pFile [open [format output_%03d/dat/Tor.%03d.dat $i_job $replica_id] w]
set qFile [open [format output_%03d/dat/img.%03d.dat $i_job $replica_id] w]
fconfigure $pFile -buffering full -buffersize 1000000
fconfigure $qFile -buffering line

load ./Rms_ext.so  Rms_ext
load ./RmsF_ext.so Rmsf_ext

set PI 3.14159265358979
set deg2Rad  [expr $PI/180.0]
set deg2Rad2 [expr $deg2Rad*$deg2Rad / 2.0] ; # There is K/2 in harmonic potential
set oldTS -1
IndF $refID
##--------------------------------------##

# ######################################## 
# ############# Tcl Forces ############### 
# ######################################## 
proc calcforces {} {
    global aid kw pFile deg2Rad deg2Rad2 Eng_b Energy sec refList wr kw_twic KwRM RMSD RMS i_job wr_HLX1 wr_HLX2 
    global phiD atmcrd ts rang IndexList KwRM_half lable RmRef cmCoor0 condR condT rangRM replica_id  kw_HLX kw_HLXTwic
    global i1 i2 i3 i4 addForce oldTS IndF refID angleFixer DelPhi Eng qFile aid2 DeltaRmsd HelixList rang_HLX
    
    set ts [getstep]
    loadcoords atmcrd
    
    set phiD   [getdihedral $atmcrd($i1) $atmcrd($i2) $atmcrd($i3) $atmcrd($i4)]
    set DelPhi [angleFixer [expr ($phiD-$wr)]]
    
    rmsD
    set DeltaRmsd [expr $RMSD - $RmRef]
    
    if {$oldTS < $ts} {
        set oldTS $ts
        puts $pFile "$phiD" 
        if {$ts % 50000 == 0} {
          flush $pFile
        }
    } else {
        IndF $refID
    }
 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# %% Apply Boundary on Helix resid 361 to 369 Phi Psi Angle %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    foreach aid1 $HelixList {
        lassign $aid1 T1 T2 T3 T4 K1 K2 K3 K4
        set phiD1  [getdihedral $atmcrd($T1) $atmcrd($T2) $atmcrd($T3) $atmcrd($T4)]
        set phiD2  [getdihedral $atmcrd($K1) $atmcrd($K2) $atmcrd($K3) $atmcrd($K4)]
        set Dalphi1 [angleFixer [expr ($phiD1-$wr_HLX1)]]
        set Dalphi2 [angleFixer [expr ($phiD2-$wr_HLX2)]]
        foreach  "inx1 inx2 inx3 inx4" $aid1 Dalpha "$Dalphi1 $Dalphi2" {
            if {abs($Dalpha) <= $rang_HLX} {
                set Eg 0.0
            } elseif  {$Dalpha > 0.0} {
                set force_b [expr -$kw_HLX*($Dalpha-$rang_HLX)*$deg2Rad]
                set Eg [addForce $force_b $kw_HLXTwic $inx1 $inx2 $inx3 $inx4]                
            } else {
                set force_b [expr -$kw_HLX*($Dalpha+$rang_HLX)*$deg2Rad]
                set Eg [addForce $force_b $kw_HLXTwic $inx1 $inx2 $inx3 $inx4]
            }
        }
    } 
 
# ######%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###### # 
# ######%%% Biased Potential on the Helix CA atoms to get RMSD %%%###### #
# ######%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###### #
    if {$condR} {
        if {$lable=="MM"} {
            rmsDF
            set Eng [expr $KwRM_half * $DeltaRmsd * $DeltaRmsd]   
        } else {
            if {abs($DeltaRmsd) <= $rangRM} {
                set Eng  0.0
            } elseif  {$DeltaRmsd > 0.0} {
                set DeltaRmsd [expr $DeltaRmsd-$rangRM]
                rmsDF
                set Eng [expr $KwRM_half * $DeltaRmsd * $DeltaRmsd]    
            } else {
                set DeltaRmsd [expr $DeltaRmsd+$rangRM]
                rmsDF
                set Eng [expr $KwRM_half * $DeltaRmsd * $DeltaRmsd]
            }
        }
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
# %% Apply Boundary condition on the Torsion Angle %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if {abs($DelPhi) <= $rang} {
            set Eng_b   0.0 
        } elseif  {$DelPhi > 0.0} {
            set force_b [expr -$kw*($DelPhi-$rang)*$deg2Rad]
            set Eng_b [addForce $force_b $kw_twic $i1 $i2 $i3 $i4] 
        } else {
            set force_b [expr -$kw*($DelPhi+$rang)*$deg2Rad]
            set Eng_b [addForce $force_b $kw_twic $i1 $i2 $i3 $i4] 
        }
# ######%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###### #
# ######%% Applying the Biased Potential on the targeted torsion angles %%###### #
# ######%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###### #
    } elseif {$condT} {  
        if {$lable=="TT"} {
            set force [expr -$kw*$DelPhi*$deg2Rad]
            set Eng   [addForce $force $kw_twic $i1 $i2 $i3 $i4]
        } else {
            if {abs($DelPhi) <= $rang} {
                set force  0.0
                set Eng    0.0
            } elseif  {$DelPhi > 0.0} {
                set force [expr -$kw*($DelPhi-$rang)*$deg2Rad]
                set Eng   [addForce $force $kw_twic $i1 $i2 $i3 $i4]
            } else {
                set force [expr -$kw*($DelPhi+$rang)*$deg2Rad]
                set Eng   [addForce $force $kw_twic $i1 $i2 $i3 $i4]
            }
        }
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %% Apply Boundary condition on the Center of mass %%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if {abs($DeltaRmsd) <= $rangRM} {
            set Eng_b   0.0 
        } elseif  {$DeltaRmsd > 0.0} {
            set DeltaRmsd [expr $DeltaRmsd-$rangRM]
            rmsDF
            set Eng_b [expr $KwRM_half * $DeltaRmsd * $DeltaRmsd] 
        } else {
            set DeltaRmsd [expr $DeltaRmsd+$rangRM]
            rmsDF
            set Eng_b [expr $KwRM_half * $DeltaRmsd * $DeltaRmsd] 
        }
    }
    
    set Energy [expr $Eng+$Eng_b]
    addenergy  $Energy
    return 
}
