
/* This file has the object container definitions */
/* Defs  **
 function Oligo()  // call new Oligo() to initialize this object
 function DoOligoCalc(form, theSequence) // direct interface with FORM
 function GetOligoMods(form) // direct interface with FORM
 function SetOligoMods(form) // direct interface with FORM
 function DoChangeInMWConcAndODOutput(form) // direct interface with FORM
 function DoOligoOutput(form, theSequenceObj) // direct interface with FORM
 function OligoCalcMWConcAndODs() // split out to make changes in famra, OD, etc easier to calc
 function DoChangeInMWConcAndODs(form) // direct interface with FORM
 function OligoCalculate()
 function OligoCount()
 function A260(choice)
 function MW(choice)
 function GC(choice)
 function Tm(choice)
 function WAKTm(choice)
 function DeltaG(choice)
 function DeltaH(choice)
 function DeltaS(choice)
 function NearestNeighborTM( choice)


 *** */

function Oligo()
{
    // initialize
    this.aCount = 0;
    this.cCount = 0;
    this.gCount = 0;
    this.tCount = 0;
    this.uCount = 0;
    this.iCount = 0;
    this.mCount = 0;
    this.rCount = 0;
    this.wCount = 0;
    this.sCount = 0;
    this.yCount = 0;
    this.kCount = 0;
    this.vCount = 0;
    this.hCount = 0;
    this.dCount = 0;
    this.bCount = 0;
    this.nCount = 0;

    this.gcValmin = 0;
    this.gcValmax = 0;
    this.mwValmin = 0;
    this.mwValmax = 0;

    this.concValmin = 0;
    this.concValmax = 0;
    this.microgramValmin = 0;
    this.microgramValmax = 0;

    this.deltaHValmin = 0;
    this.deltaGValmin = 0;
    this.deltaSValmin = 0;
    this.deltaHValmax = 0;
    this.deltaGValmax = 0;
    this.deltaSValmax = 0;

    this.RlogK = 0;

    this.basicTmValmin = 0;
    this.adjustedTmValmin = 0;
    this.nearestNeighborTmValmin = 0;
    this.basicTmValmax = 0;
    this.adjustedTmValmax = 0;
    this.nearestNeighborTmValmax = 0;

    this.eA_A260min=0;
    this.eC_A260min=0;
    this.eG_A260min=0;
    this.eT_A260min=0;
    this.eU_A260min=0;

    this.eA_A260max=0;
    this.eC_A260max=0;
    this.eG_A260max=0;
    this.eT_A260max=0;
    this.eU_A260max=0;

    //effective counts for MW()
    this.eA_MWmin=0;
    this.eC_MWmin=0;
    this.eG_MWmin=0;
    this.eT_MWmin=0;
    this.eU_MWmin=0;
    this.eI_MWmin=0;

    this.eA_MWmax=0;
    this.eC_MWmax=0;
    this.eG_MWmax=0;
    this.eT_MWmax=0;
    this.eU_MWmax=0;
    this.eI_MWmax=0;
    this.eGC_min=0;
    this.eGC_max=0;


    this.aaCount = 0;
    this.atCount = 0;
    this.taCount = 0;
    this.caCount = 0;
    this.gtCount = 0;
    this.ctCount = 0;
    this.gaCount = 0;
    this.cgCount = 0;
    this.gcCount = 0;
    this.ggCount = 0;
    this.IUpairVals_min=new Array(0,0,0);
    this.IUpairVals_max=new Array(0,0,0);
    this.hasIUpacBase=0;
    this.isDeoxy=1;
    this.isSingleStranded=1;

    // arrays
    this.seqArray;
    this.revSeqArray;

    // function associations
    this.Tm = Tm;
    this.WAKTm = WAKTm;
    this.GC = GC;
    this.MW = MW;
    this.A260 = A260;
    this.DeltaH = DeltaH;
    this.DeltaG = DeltaG;
    this.DeltaS = DeltaS;
    this.NearestNeighborTM = NearestNeighborTM;
    this.OligoCount = OligoCount;
    this.OligoCalculate = OligoCalculate;
    this.OligoCalcMWConcAndODs=OligoCalcMWConcAndODs;
    this.DoChangeInMWConcAndODs = DoChangeInMWConcAndODs; //public interface
    this.DoChangeInMWConcAndODOutput = DoChangeInMWConcAndODOutput; //public interface
    this.DoOligoCalc = DoOligoCalc;  //public interface
    this.GetOligoMods = GetOligoMods;  //public interface
    this.SetOligoMods = SetOligoMods;  //public interface
    this.DoOligoOutput = DoOligoOutput;  //public interface

    // publically set values
    this.Sequence="";
    this.revSequence="";
    this.ODs = 1;
    this.FivePrimeModification = "";
    this.FivePrimeMW = 0;
    this.ThreePrimeModification = "";
    this.ThreePrimeMW = 0;

    //modifications
    this.ThreePrimeA260mod = 0;
    this.FivePrimeA260mod = 0;
    //this.hexCount = 0;
    //this.tamraCount = 0;
    this.saltConcentration = 50;
    this.primerConcentration = 50;
}

// direct interface with FORM

function GetOligoMods(form)  {
    this.FivePrimeModification = form.FivePrime.options[form.FivePrime.selectedIndex].text;
    this.FivePrimeMW = parseInt( form.FivePrime.options[form.FivePrime.selectedIndex].value);
    this.ThreePrimeModification = form.ThreePrime.options[form.ThreePrime.selectedIndex].text;
    this.ThreePrimeMW = parseInt(form.ThreePrime.options[form.ThreePrime.selectedIndex].value);
    /*
     this.famCount   = form.fam.value;
     this.hexCount   = form.hex.value;
     this.tetCount   = form.tet.value;
     this.tamraCount = form.tamra.value;
     */
}

// direct interface with FORM

function SetOligoMods(form)  {
    form.FivePrime.options.selectedIndex=0;
    for (var i = 0; i < form.FivePrime.options.length; i++) {
        if (this.FivePrimeModification.length > 0) {
            if (form.FivePrime.options[i].text.indexOf(this.FivePrimeModification) >= 0) {
                form.FivePrime.options[i].selected=true;
                form.FivePrime.options.selectedIndex=i;
            } else {
                form.FivePrime.options[i].selected=false;
            }
        }
    }
    form.ThreePrime.options.selectedIndex=0;
    for (var i = 0; i < form.ThreePrime.options.length; i++) {
        if (this.ThreePrimeModification.length > 0) {
            if (form.ThreePrime.options[i].text.indexOf(this.ThreePrimeModification) >= 0) {
                form.ThreePrime.options[i].selected=true;
                form.ThreePrime.options.selectedIndex=i;
            } else {
                form.ThreePrime.options[i].selected=false;
            }
        }
    }
}

// direct interface with FORM

function DoOligoCalc(form, theSequence)  {
    temp=CheckBase(theSequence); // external call
    if (temp==-1) {return;} // do not continue if the sequence is not valid!
    this.Sequence=temp;
    this.ODs        = form.ODs.value;
    this.saltConcentration = parseFloat(form.saltConcBox.value); // parseInt(form.saltConcBox.value);
    this.primerConcentration = parseFloat(form.primerConcBox.value);  // parseInt(form.primerConcBox.value);
    this.hasIUpacBase=AreThereIUpacBases(this.Sequence);
    if (debug) alert("Sequence length="+this.Sequence.length+"; Sequence="+this.Sequence);
    temp = form.deoxy.options[form.deoxy.selectedIndex].value;
    if (temp.indexOf("DNA") >= 0) {
        this.isDeoxy	= 1;
    } else {
        this.isDeoxy	= 0;
    }
    if (temp.indexOf("ds") >= 0) {
        this.isSingleStranded	= 0;
    } else {
        this.isSingleStranded	= 1;
    }
    this.OligoCalculate();
}

// direct interface with FORM
function DoChangeInMWConcAndODOutput(form)  {
    if (debug) alert('MW='+this.mwValmin);
    form.ODreadOnly.value = this.ODs;
    if(!this.hasIUpacBase) {
        if (this.isSingleStranded) {
            form.mwBox.value = this.mwValmin;
        } else {
            form.mwBox.value=theOligo.mwValmin+theComplement.mwValmin;
        }

        form.micromolarConc.value = this.concValmin;
        form.micrograms.value = this.microgramValmin;
    } else {
        if (this.isSingleStranded) {
            form.mwBox.value = this.mwValmin+" to "+this.mwValmax;
        } else {
            var tempMIn=theOligo.mwValmin+theComplement.mwValmin;
            var tempMax=theOligo.mwValmax+theComplement.mwValmax;
            form.mwBox.value = tempMIn+" to "+tempMax;
        }

        form.micromolarConc.value=this.concValmin+" to "+this.concValmax;
        form.micrograms.value = this.microgramValmin+" to "+this.microgramValmax;
    }
}

// direct interface with FORM
function DoOligoOutput(form, theSequenceObj)  {
    this.DoChangeInMWConcAndODOutput(form);
    if(!this.hasIUpacBase){
        form.gcBox.value    = this.gcValmin;
        form.tmBox.value    = this.basicTmValmin;
        form.WAKtmBox.value = this.adjustedTmValmin;
        form.nTmBox.value   = this.nearestNeighborTmValmin;
    }else{
        form.gcBox.value    = this.gcValmin+" to "+this.gcValmax;
        form.tmBox.value    = this.basicTmValmin+" to "+this.basicTmValmax;
        form.WAKtmBox.value = this.adjustedTmValmin+" to "+this.adjustedTmValmax;
        form.nTmBox.value   = this.nearestNeighborTmValmin+" to "+this.nearestNeighborTmValmax;
    }
    form.RlogKBox.value =this.RlogK;
    if(!this.hasIUpacBase) {
        form.deltaHBox.value = this.deltaHValmin;
        form.deltaGBox.value = this.deltaGValmin;
        form.deltaSBox.value = this.deltaSValmin;
    } else {
        form.deltaHBox.value = this.deltaHValmin+" to "+this.deltaHValmax;
        form.deltaGBox.value = this.deltaGValmin+" to "+this.deltaGValmax;
        form.deltaSBox.value = this.deltaSValmin+" to "+this.deltaSValmax;
    }
    form.lBox.value = this.Sequence.length;
    theSequenceObj.value =FormatBaseString(this.Sequence);
}

// split out to make changes in famra, OD, etc easier to calc
function OligoCalcMWConcAndODs()  {
    /*** Now do MW calculation ***/
    if(!this.hasIUpacBase) {
        this.mwValmin = this.MW("min");
        this.mwValmax = this.mwValmin;
    } else {
        this.mwValmin = this.MW("min");
        this.mwValmax = this.MW("max");
    }
    if (debug) alert("oligoCalc mw val="+this.mwValmin);
    /*** Now do A260 concentration calculation ***/
    this.concValmin = this.A260("min");
    this.concValmax = this.A260("max");
    this.microgramValmin = micrograms(this.mwValmin, this.concValmin); // external call
    this.microgramValmax = micrograms(this.mwValmax, this.concValmax);
}

// direct interface with FORM
function DoChangeInMWConcAndODs(form)  {
    this.FivePrimeModification = form.FivePrime.options[form.FivePrime.selectedIndex].text;
    this.FivePrimeMW = parseInt(form.FivePrime.options[form.FivePrime.selectedIndex].value);
    this.ThreePrimeModification = form.ThreePrime.options[form.ThreePrime.selectedIndex].text;
    this.ThreePrimeMW = parseInt(form.ThreePrime.options[form.ThreePrime.selectedIndex].value);
    this.FivePrimeA260mod = ParseA260modification(this.FivePrimeModification);
    this.ThreePrimeA260mod = ParseA260modification(this.ThreePrimeModification);
    /*
     this.famCount   = form.fam.value;
     this.hexCount   = form.hex.value;
     this.tetCount   = form.tet.value;
     this.tamraCount = form.tamra.value;
     */
    this.ODs        = form.ODs.value;
    temp = form.deoxy.options[form.deoxy.selectedIndex].value;
    if (temp.indexOf("DNA") >= 0) {
        this.isDeoxy	= 1;
    } else {
        this.isDeoxy	= 0;
    }
    if (temp.indexOf("ds") >= 0) {
        this.isSingleStranded	= 0;
    } else {
        this.isSingleStranded	= 1;
    }
    this.OligoCalculate();
}


function OligoCalculate() {
    this.OligoCount();
    /*** Do GC calculation ***/
    if(!this.hasIUpacBase) {
        this.gcValmin = this.GC("min");
        this.gcValmax = this.gcValmin;
    } else {
        this.gcValmin = this.GC("min");
        this.gcValmax = this.GC("max");
    }
    this.OligoCalcMWConcAndODs()

    /** Calculate numbers for Thermodynamic TM calculation **/
    if(!this.hasIUpacBase)
    {
        this.deltaHValmin = this.DeltaH("min");
        this.deltaGValmin = this.DeltaG("min");
        this.deltaSValmin = this.DeltaS("min");
        this.deltaHValmax = this.deltaHValmin;
        this.deltaGValmax = this.deltaGValmin;
        this.deltaSValmax = this.deltaSValmin;
    } else {
        this.deltaHValmin = this.DeltaH("min");
        this.deltaGValmin = this.DeltaG("min");
        this.deltaSValmin = this.DeltaS("min");
        this.deltaHValmax = this.DeltaH("max");
        this.deltaGValmax = this.DeltaG("max");
        this.deltaSValmax = this.DeltaS("max");
    }
    if(!this.hasIUpacBase){
        this.basicTmValmin           = this.Tm("min");
        this.adjustedTmValmin        = this.WAKTm("min");
        this.nearestNeighborTmValmin = this.NearestNeighborTM("min");
        this.basicTmValmax           = this.basicTmValmin;
        this.adjustedTmValmax        = this.adjustedTmValmin;
        this.nearestNeighborTmValmax = this.adjustedTmValmin;
    }else{
        this.basicTmValmin           = this.Tm("min");
        this.adjustedTmValmin        = this.WAKTm("min");
        this.nearestNeighborTmValmin = this.NearestNeighborTM("min");
        this.basicTmValmax           = this.Tm("max");
        this.adjustedTmValmax        = this.WAKTm("max");
        this.nearestNeighborTmValmax = this.NearestNeighborTM("max");
    }
}

function OligoCount(){
    this.aCount = CountChar("A",this.Sequence);
    this.cCount = CountChar("C",this.Sequence);
    this.gCount = CountChar("G",this.Sequence);
    this.tCount = CountChar("T",this.Sequence);
    this.uCount = CountChar("U",this.Sequence);
    this.iCount = CountChar("I",this.Sequence);
    this.mCount = CountChar("M",this.Sequence);
    this.rCount = CountChar("R",this.Sequence);
    this.wCount = CountChar("W",this.Sequence);
    this.sCount = CountChar("S",this.Sequence);
    this.yCount = CountChar("Y",this.Sequence);
    this.kCount = CountChar("K",this.Sequence);
    this.vCount = CountChar("V",this.Sequence);
    this.hCount = CountChar("H",this.Sequence);
    this.dCount = CountChar("D",this.Sequence);
    this.bCount = CountChar("B",this.Sequence);
    this.nCount = CountChar("N",this.Sequence);
    //Effective a,c,g,t count for different calculations
    //effective counts for A260min()
    this.eA_A260min=this.aCount;
    this.eU_A260min=this.uCount;
    this.eI_A260min=this.iCount;
    this.eC_A260min=this.cCount+this.mCount+this.sCount+this.yCount+this.vCount
        +this.hCount+this.bCount+this.nCount;
    this.eG_A260min=this.gCount+this.rCount;
    this.eT_A260min=this.tCount+this.wCount+this.kCount+this.dCount;

    this.eA_A260max=this.aCount+this.mCount+this.rCount+this.wCount+this.vCount
        +this.hCount+this.bCount+this.nCount;
    this.eU_A260max=this.uCount;
    this.eC_A260max=this.cCount;
    this.eG_A260max=this.gCount+this.sCount+this.kCount+this.dCount;
    this.eT_A260max=this.tCount+this.yCount;

    //effective counts for MW()
    this.eA_MWmin=this.aCount+this.rCount;
    this.eC_MWmin=this.eC_A260min;
    this.eU_MWmin=this.eU_A260min;
    this.eI_MWmin=this.eI_A260min;
    this.eG_MWmin=this.gCount;
    this.eT_MWmin=this.eT_A260min;

    this.eA_MWmax=this.aCount+this.mCount+this.wCount+this.hCount;
    this.eC_MWmax=this.cCount;
    this.eG_MWmax=this.gCount+this.rCount+this.sCount+this.kCount+this.vCount+
        this.dCount+this.bCount+this.nCount;
    this.eT_MWmax=this.tCount+this.yCount;

    this.eGC_min=this.gCount + this.cCount +this.sCount;
    this.eGC_max=(this.Sequence.length)-this.aCount-this.tCount-this.wCount-this.uCount;

    // count Nearest Neighbors
    this.aaCount = CountNeighbors("AA",this.Sequence)+CountNeighbors("TT",this.Sequence)+CountNeighbors("UU",this.Sequence);
    this.atCount = CountNeighbors("AT",this.Sequence)+CountNeighbors("AU",this.Sequence);
    this.taCount = CountNeighbors("TA",this.Sequence)+CountNeighbors("UA",this.Sequence);
    this.caCount = CountNeighbors("CA",this.Sequence)+CountNeighbors("TG",this.Sequence)+CountNeighbors("UG",this.Sequence);
    this.gtCount = CountNeighbors("GT",this.Sequence)+CountNeighbors("AC",this.Sequence)+CountNeighbors("GU",this.Sequence);
    this.ctCount = CountNeighbors("CT",this.Sequence)+CountNeighbors("AG",this.Sequence)+CountNeighbors("CU",this.Sequence);
    this.gaCount = CountNeighbors("GA",this.Sequence)+CountNeighbors("TC",this.Sequence)+CountNeighbors("UC",this.Sequence);
    this.cgCount = CountNeighbors("CG",this.Sequence);
    this.gcCount = CountNeighbors("GC",this.Sequence);
    this.ggCount = CountNeighbors("GG",this.Sequence)+CountNeighbors("CC",this.Sequence);
    //Calculate IUpac pairs
    /*----08/02/99 fix one bug for nearest Neighbor Calc, */
    for(var j=0; j<3; j++)
    {
        this.IUpairVals_min[j]=0;
        this.IUpairVals_max[j]=0;
    }
    /*******************************************************/
    for(var i=1; i<this.Sequence.length;i++) {	//first base can not be IUpacbase
        var base0=this.Sequence.charAt((i-1));
        var base=this.Sequence.charAt(i);
        var temp=new Array(0,0,0);

        temp=CalcIUpair(base0, base, i, this.Sequence, "min");
        if (debug) alert("temp "+temp[0]+" "+temp[1]+" "+temp[2]);
        for(var j=0; j<3; j++) {
            this.IUpairVals_min[j]+=temp[j];
        }
        if (debug) alert("mim"+this.IUpairVals_min[0]+" "+this.IUpairVals_min[1]+" "+this.IUpairVals_min[2]);
        temp=CalcIUpair(base0, base, i, this.Sequence, "max");
        if (debug) alert("temp "+temp[0]+" "+temp[1]+" "+temp[2]);
        for( j=0; j<3; j++) {
            this.IUpairVals_max[j]+=temp[j];
        }
        if (debug) alert("max"+this.IUpairVals_max[0]+" "+this.IUpairVals_max[1]+" "+this.IUpairVals_max[2]);
    }
}

function ParseA260modification(text) {
    if ( /Amino dT/.test(text) )
        return 8700;
    if ( /^BHQ/.test(text) )
        return 8000;
    if ( /^Carboxy-dT/.test(text) )
        return 14700;
    if ( /^Cyanine (CY) 3.5/.test(text) )
        return 24000;
    if ( /^Cyanine (CY) 5.5/.test(text) )
        return 21500;
    if ( /^Cyanine (CY) 3/.test(text) )
        return 4900;
    if ( /^Cyanine (CY) 5/.test(text) )
        return 10000;
    if ( /^Dabsyl/.test(text) )
        return 32000;
    if ( /6-FAM/.test(text) )
        return 20960;
    if ( /Hexachloro-Fluorescein/.test(text) )
        return 31580;
    if ( /^TET/.test(text) )
        return 16255;
    if ( /Fluorescein/.test(text) )
        return 13700;
    if ( /^Inosine/.test(text) )
        return 7200;
    if ( /dCytidine/.test(text) )
        return 3300;
    if ( /^JOE/.test(text) )
        return 75000;
    if ( /^Oregon Green/.test(text) )
        return 76000;
    if ( /^Rhodamine Fluorophore/.test(text) )
        return 92500;
    if ( /^Rhodamine-Green/.test(text) )
        return 68000;
    if ( /^6-Carboxy-X-rhodamine/.test(text) )
        return 82000;
    if ( /TAMRA/.test(text) )
        return 31980;
    if ( /^Texas Red/.test(text) )
        return 80000;

    return 0;
}

/* Molar Absorptivity
 Values from ABI, units M-1 cm-1
 A	15200	G	12010	C	7050	T	8400	6' FAM	20960	TET	16255
 HEX	31580	TAMRA	31980
 Need ribonucleotide molar absorptivities
 U 20800
 Assume 1 OD of a standard 1ml solution (1 cm pathlength)
 -----  */
function A260(choice) {
    var div;
    if (this.Sequence.length > 0) {
        if (this.isSingleStranded) {
            if (this.isDeoxy) {
                if(choice=="min"){	//calculates the minimum value, use the e_A260max since it is the divident
                    div=(this.eA_A260max * 15200 +
                        this.eG_A260max * 12010 +
                        this.eC_A260max *  7050 +
                        this.eT_A260max *  8400 +
                        this.eU_A260max * 9800 +
                        this.ThreePrimeA260mod +
                        this.FivePrimeA260mod);
                }else{ //choice=="max"
                    div=(this.eA_A260min * 15200 +
                        this.eG_A260min * 12010 +
                        this.eC_A260min *  7050 +
                        this.eT_A260min *  8400 +
                        this.eU_A260min * 9800 +
                        this.ThreePrimeA260mod +
                        this.FivePrimeA260mod);
                }
            } else {
                if(choice=="min"){	//calculates the minimum value, use the e_A260max since it is the divident
                    div=(this.eA_A260max * 15400 +
                        this.eG_A260max * 13700 +
                        this.eC_A260max *  9000 +
                        this.eT_A260max *  9400 +
                        this.eU_A260max * 10000 +
                        this.ThreePrimeA260mod +
                        this.FivePrimeA260mod);
                }else{ //choice=="max"
                    div=(this.eA_A260min * 15400 +
                        this.eG_A260min * 13700 +
                        this.eC_A260min *  9000 +
                        this.eT_A260min *  9400 +
                        this.eU_A260min * 10000 +
                        this.ThreePrimeA260mod +
                        this.FivePrimeA260mod);
                }
            }
        } else {
            if (this.isDeoxy) {
                div=this.Sequence.length * 2 * 6400;
            } else {
                div=this.Sequence.length * 2 * 8000;
            }
        }
        // units are in microMoles/liter
        return (Math.round( this.ODs * 	1000000000/div)/1000);

    }
    return "";
}

function GetMod_A260() {
    return 0;
}

function MW(choice) {
    var mw;
    if (this.Sequence.length > 0) {
        if (this.isDeoxy) {
            if(choice=="min"){
                mw=     313.21 * this.eA_MWmin +
                    329.21 * this.eG_MWmin +
                    289.18 * this.eC_MWmin +
                    304.2  * this.eT_MWmin +
                    290.169 * this.eU_MWmin +
                    314. * this.eI_MWmin
                    - 61.96;
            }else{
                mw=     313.21 * this.eA_MWmax +
                    329.21 * this.eG_MWmax +
                    289.18 * this.eC_MWmax +
                    304.2 * this.eT_MWmax +
                    290.169 * this.eU_MWmax +
                    314. * this.eI_MWmax
                    - 61.96;
            }
        } else { // is riboNucleotide
            if(choice=="min"){
                mw=     329.21 * this.eA_MWmin +
                    345.21 * this.eG_MWmin +
                    305.18 * this.eC_MWmin +
                    320.2  * this.eT_MWmin +
                    306.169 * this.eU_MWmin +
                    330  * this.eI_MWmin
                    +159;

            }else{
                mw=     329.21 * this.eA_MWmax +
                    345.21 * this.eG_MWmax +
                    305.18 * this.eC_MWmax +
                    320.2 * this.eT_MWmax +
                    306.169 * this.eU_MWmax +
                    330  * this.eI_MWmax
                    +159;
            }
        }
        var new_mw = mw+this.FivePrimeMW;
        new_mw += this.ThreePrimeMW;
        new_mw= Math.round( 10 * new_mw)/10;
        return new_mw;
    }
    return "";
}

function GC(choice)
{
    if (this.Sequence.length > 0) {
        if(choice=="min"){
            return Math.round( 100 * this.eGC_min /this.Sequence.length )
        }else{
            return Math.round(100*this.eGC_max/this.Sequence.length);
        }
    }
    return "";
}

function Tm(choice)
{
    if (this.Sequence.length > 0) {
        if (this.Sequence.length < 14) {
            if(choice=="min"){
                return Math.round((2 * (this.Sequence.length-this.eGC_min) + 4 * (this.eGC_min))*10)/10;
            }else{
                return Math.round((2 * (this.Sequence.length-this.eGC_max) + 4 * (this.eGC_max))*10)/10;
            }
        }
        else {
            if(choice=="min"){
                return Math.round((64.9 +  41*((this.eGC_min - 16.4) / this.Sequence.length))*10)/10;
            }else{
                return Math.round((64.9 + 41*((this.eGC_max - 16.4) / this.Sequence.length))*10)/10;
            }
        }
    }
    return "";
}

function WAKTm(choice){
    if (this.Sequence.length > 0) {
        if (this.isDeoxy) {
            if (this.Sequence.length < 14) {
                if(choice=="min"){

                    return  Math.round((2 * (this.Sequence.length-this.eGC_min) + 4 * (this.eGC_min)+
                        21.6+(7.21*Math.log(this.saltConcentration/1000)))*10)/10;
                } else{
                    return  Math.round((2 * (this.Sequence.length-this.eGC_max) + 4 * (this.eGC_max)+
                        21.6+(7.21*Math.log(this.saltConcentration/1000)))*10)/10;
                }
            }
            else {
                if (choice=="min"){
                    return Math.round((100.5 + (0.41*this.gcValmin) - (820 / this.Sequence.length)+
                        (7.21*Math.log(this.saltConcentration/1000)))*10)/10;
                } else {
                    return Math.round((100.5 + (0.41*this.gcValmax) - (820 / this.Sequence.length)+
                        (7.21*Math.log(this.saltConcentration/1000)))*10)/10;
                }
            }
        } else {
            if (this.Sequence.length > 20) {
                if (choice=="min"){
                    return Math.round((79.8 + (0.584*this.gcValmin) + (11.8*(this.gcValmin/100)*(this.gcValmin/100)) - (820 / this.Sequence.length)+
                        (8.03*Math.log(this.saltConcentration/1000)))*10)/10;
                } else {
                    return Math.round((79.8 + (0.584*this.gcValmax) - (820 / this.Sequence.length)+
                        (8.03*Math.log(this.saltConcentration/1000)))*10)/10;
                }
            } else {
                return "too short";
            }
        }
    }
    return "";
}

/* ******
 deltaH=deltaG + TdeltaS
 delta G =  RT ln ([DNA-primer complex]/([DNA][primer]))
 T = deltaH/(deltaS + R ln([DNA-primer complex]/([DNA][primer])))  + 7.21* ln [Na concentration])

 Assume [DNA-primer complex] = [DNA]

 NOTE: Javascript Math.log() calls a NATURAL LOG, not Log base 10!
 */

function DeltaG(choice)
{
    if (this.Sequence.length > 7) {
        var val= -5.0;
        // Helix initiation Free Energy of 5 kcal.
        // symmetry function: if symmetrical, subtract another 0.4
        val+=1.2*this.aaCount;
        val+=0.9*this.atCount;
        val+=0.9*this.taCount;
        val+=1.7*this.caCount;
        val+=1.5*this.gtCount;
        val+=1.5*this.ctCount;
        val+=1.5*this.gaCount;
        val+=2.8*this.cgCount;
        val+=2.3*this.gcCount;
        val+=2.1*this.ggCount;
        if(choice=="min"){
            val+=this.IUpairVals_min[0];
        }else{
            val+=this.IUpairVals_max[0];
        }
        return Math.round((1000*val))/1000;
    }
    return "";
}

function DeltaH(choice)
{
    if (this.Sequence.length > 7) {
        var val= 0.0;
        if (this.isDeoxy) {
            val+=8.0*this.aaCount;
            val+=5.6*this.atCount;
            val+=6.6*this.taCount;
            val+=8.2*this.caCount;
            val+=9.4*this.gtCount;
            val+=6.6*this.ctCount;
            val+=8.8*this.gaCount;
            val+=11.8*this.cgCount;
            val+=10.5*this.gcCount;
            val+=10.9*this.ggCount;
        } else {
            val+=6.8*this.aaCount;
            val+=9.38*this.atCount;
            val+=7.69*this.taCount;
            val+=10.44*this.caCount;
            val+=11.4*this.gtCount;
            val+=10.48*this.ctCount;
            val+=12.44*this.gaCount;
            val+=10.64*this.cgCount;
            val+=14.88*this.gcCount;
            val+=13.39*this.ggCount;
        }
        if(choice=="min"){
            val+=this.IUpairVals_min[1];
        }else{
            val+=this.IUpairVals_max[1];
        }
        return Math.round((1000*val))/1000;
    }
    return "";
}

/*
 RNA nearest neighbor information

 Xia T., SantaLucia J., Burkard M.E., Kierzek R., Schroeder S.J., Jiao X.,
 Cox C., Turner D.H. (1998). Thermodynamics parameters for an expanded
 nearest-neighbor model for formation of RNA duplexes with Watson-Crick
 base pairs. Biochemistry 37: 14719-14735

 DeltaH   DeltaS
 J.mol-1  J.mol-1.K-1

 AA  -6820.0  -19.0
 AC -10200.0  -26.2
 AG  -7600.0  -19.2
 AU  -9380.0  -26.7
 CA -10440.0  -26.9
 CC -12200.0  -29.7
 CG -10640.0  -26.7
 CU -10480.0  -27.1
 GA -12440.0  -32.5
 GC -14880.0  -36.9
 GG -13390.0  -32.7
 GU -11400.0  -29.5
 UA  -7690.0  -20.5
 UC -13300.0  -35.5
 UG -10500.0  -27.8
 UU  -6600.0  -18.4
 IA   3720.0  -10.5
 IG   3610.0    1.5
 */


function DeltaS(choice)
{
    if (this.Sequence.length > 7) {
        var val=0;
        if (this.isDeoxy) {
            val+=21.9*this.aaCount;
            val+=15.2*this.atCount;
            val+=18.4*this.taCount;
            val+=21.0*this.caCount;
            val+=25.5*this.gtCount;
            val+=16.4*this.ctCount;
            val+=23.5*this.gaCount;
            val+=29.0*this.cgCount;
            val+=26.4*this.gcCount;
            val+=28.4*this.ggCount;
        } else {
            val+=19.0*this.aaCount;
            val+=26.7*this.atCount;
            val+=20.5*this.taCount;
            val+=26.9*this.caCount;
            val+=29.5*this.gtCount;
            val+=27.1*this.ctCount;
            val+=32.5*this.gaCount;
            val+=26.7*this.cgCount;
            val+=36.9*this.gcCount;
            val+=32.7*this.ggCount;
        }
        if(choice=="min"){
            val+=this.IUpairVals_min[2];
        }else{
            val+=this.IUpairVals_max[2];
        }
        return Math.round((1000*val))/1000;
    }
    return "";
}

function NearestNeighborTM( choice) {
    var theReturn = "";
    if (this.Sequence.length > 7) {
        //
        var K = 1/(this.primerConcentration*1E-9);  // Convert from nanomoles to moles
        var R = 1.987;  //cal/(mole*K);
        var RlnK = R*Math.log(K); // javascript log is the natural log (ln)
        this.RlogK =Math.round(1000* RlnK)/1000;
        // Helix initiation Free Energy of 3.4 kcal (Sugimoto et al, 1996)
        // symmetry function: if symmetrical, subtract another 0.4 - not implemented
        if (choice =="min") {
            theReturn = 1000*((this.DeltaH("min")-3.4)/(this.DeltaS("min")+RlnK));
            theReturn += -272.9;  // Kelvin to Centigrade
            theReturn += 7.21*Math.log (this.saltConcentration/1000);
            theReturn = Math.round(100*theReturn)/100;
        } else {
            theReturn = 1000*((this.DeltaH("max")-3.4)/(this.DeltaS("max")+RlnK));
            theReturn += -272.9; // Kelvin to Centigrade
            theReturn += 7.21*Math.log (this.saltConcentration/1000);
            theReturn = Math.round(100*theReturn)/100;
        }
    } else {
        this.RlogK ="";
    }
    return theReturn;
}

