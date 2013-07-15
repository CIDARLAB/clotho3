/**
 * Global Variables
 **/

var calculators = new Array();
var timer;
var mostRecentDnaproteinField;
var mostRecentTempField;
var isInitialRequest = true;


/**
 * Document.ready()
 **/
$(function() {

    //Test for IE7 or IE8, display only unordered list, never select box
    if (navigator.appVersion.indexOf("MSIE 7.") != -1 || navigator.appVersion.indexOf("MSIE 8.") != -1) {
        $('#calculatorList').removeClass('visible-tablet').removeClass('visible-desktop');
        $('#calculatorSelectDiv').removeClass('visible-phone');
        $('#calculatorSelect').removeClass('visible-phone');
    }

    if (!('ontouchstart' in document.documentElement)) {
        $('li.listCalc').add('li.topListItem').add('li.bottomListItem').hover(function(){
            $(this).css("color", "rgb(253,183,19)");
        }, function(){
            $(this).css("color", "#ffffff");
        });
    }  else {
        $('li.listCalc').add('li.topListItem').add('li.bottomListItem').bind('touchstart', function(){
            $(this).css("color", "rgb(253,183,19)");
        });
        
        $('li.listCalc').add('li.topListItem').add('li.bottomListItem').bind('touchend', function(){
            $(this).css("color", "#ffffff");    
        });
    }
    createCalculatorObjects();
    
    //Bind handler to select menu for phone-sized window
    $('#calculatorSelect').change(function() {
        var calcNumber = $(this).find(":selected").val();
        
        if (calcNumber === '0') {
            $('#centerPanel').empty();
            $('#rightPanel').empty();
        } else {
            calculators[calcNumber-1].loadCalculator();
        }
    });
    
    var rerenderFormula = function() {
        var formulaDiv = document.getElementById('formulaDiv');
        MathJax.Hub.Queue(["Rerender", MathJax.Hub, formulaDiv]); 
    }

    $(window).resize(function() {
        clearTimeout(timer);
        timer = setTimeout(rerenderFormula, 500);
    });
    //Load first calculator
    var param = getParameterByName('calc');

    if (param ==='') {
        calculators[0].loadCalculator();
    } else {
        var calculator;
        for (var i = 0; i < calculators.length; i++) {
            if(calculators[i].id === param) {
                calculator = calculators[i];
            }
        };
        calculator.loadCalculator();
    }
});

// getParameterByName
function getParameterByName(name) {
    name = name.replace(/[\[]/, "\\\[").replace(/[\]]/, "\\\]");
    var regexS = "[\\?&]" + name + "=([^&#]*)";
    var regex = new RegExp(regexS);
    var results = regex.exec(window.location.href);
    if ( results == null ) {
        return "";
    } else {
        return decodeURIComponent(results[1].replace(/\+/g, " "));
    }
}

/**
 * Calculator Object Creation, Contsructor &  Prototype Methods
 **/
function clearAnswerWarningField(e) {
    var parentID = "#"+$(e.target).closest(".calculator").attr('id');
    $('.calcField').each(function(){
        $(this).removeClass('error info warning');
        $(this).children('.controls').children('.help-inline').html("");
    });
    $(parentID+"Answer").hide().html("").removeClass('answer');
    $('.messageHolder').hide();
}


function createCalculatorObjects()
{
    var calculatorIdArray = new Array("ugpmols", "pmolng", "ugmlpmolul", "pmolulugml", "pmolends", "ratio", "odConvert", "kdapmolug", "dnaprotein", "dilution", "molarity", "temp", "tm");
    var calculatorNameArray = new Array("dsDNA: &micro;g to pmol",
                                        "dsDNA: pmol to &micro;g",
                                        "ssDNA: Micrograms per Milliliter to Picomoles per Microliter",
                                        "ssDNA: Picomoles per Microliter to Micrograms per Milliliter",
                                        "Linear DNA: &micro;g to pmol of Ends",
                                        "Ligations: Molar Ratio of Insert:Vector",
                                        "Nucleic Acid: OD<sub>260</sub> to &micro;g/ml",
                                        "Molar Conversion",
                                        "Coding Capacity of DNA",
                                        "Dilution",
                                        "Molarity",
                                        "Temperature Conversion",
                                        "T<sub>m</sub> for Oligos");
    
    //Clear function used by majority of calculators

    var clearFields =   function(e){
                            var parentID = "#"+$(e.target).closest(".calculator").attr('id');
                            clearAnswerWarningField(e);
                            if(parentID=="#tm"){
                                sequence.value="";
                                var myselect = $("select#primers");
                                myselect[0].selectedIndex = 0;
                                return;
                            }
                            $(parentID+" form input").val('');
                        }
        
    for(var i = 0; i < calculatorIdArray.length; i++)
    {
        //Obtain number of the calculator in 2-digit format - used for retrieving formula image file.
        var number = (i + 1) + '';
        
        if (number.length < 2)
        {
            number = '0' + number;
        }
        
        //Assign calculation and standard clear functions to calculator object        
        var calcFunction;
        var clearFunction = clearFields;
        
        switch(calculatorIdArray[i]) {
        case "ugpmols":
            calcFunction =  function(e){
                                clearAnswerWarningField(e);
                                var inputFields = $('form').find('input');
                                if (validateInputValues(inputFields)) {
                                    var answer = (parseFloat($('#ugDNA').val()) / parseFloat($('#n').val())) * 1515.1;
                                    answer = biomathRounding(answer);
                                    var units = "pmols of DNA";
                                    var parentID = '#ugpmols';
                                    returnAnswer(answer,parentID,units);
                                };
                            }
            break;
        case "pmolng":
            calcFunction =  function (e){
                                clearAnswerWarningField(e);
                                var inputFields = $('form').find('input');
                                if (validateInputValues(inputFields)) {
                                    var answer = (parseFloat($('#pmol').val()) * parseFloat($('#bp').val())) * 0.00066;
                                    answer = biomathRounding(answer);
                                    var units = "&micro;g of DNA";
                                    var parentID = '#pmolng';
                                    returnAnswer(answer,parentID,units);
                                }
                            }
            break;
        case "ugmlpmolul":
            calcFunction =  function (e){
                                clearAnswerWarningField(e);
                                //var inputFields = [$('#bases'), $('#conc')];
                                var inputFields = $('form').find('input');
                                if (validateInputValues(inputFields)) {
                                    var answer = parseFloat($("#conc").val()) * 3.03 / parseFloat($("#bases").val());
                                    answer = biomathRounding(answer);
                                    var units = "pmol/&micro;l of DNA";
                                    var parentID = '#ugmlpmolul';
                                    returnAnswer(answer,parentID,units);
                                }
                            }
            break;
        case "pmolulugml":
            calcFunction =  function (e){
                                clearAnswerWarningField(e);
                                var inputFields = $('form').find('input');
                                if (validateInputValues(inputFields)) {
                                    var answer = parseFloat($("#oligoConc").val()) * 0.33 * parseFloat($("#oligoLength").val());
                                    answer = biomathRounding(answer);
                                    var units = "&micro;g/ml of DNA";
                                    var parentID = '#pmolulugml';
                                    returnAnswer(answer,parentID,units);
                                };
                                
                            }
            break;
        case "pmolends":
            calcFunction =  function (e){
                                clearAnswerWarningField(e);
                                var inputFields = $('form').find('input');
                                if (validateInputValues(inputFields)) {
                                    var answer = parseFloat($("#linearUg").val()) * 3.03 / parseFloat($("#linearKb").val());
                                    answer = biomathRounding(answer);
                                    var units = "pmol ends";
                                    var parentID = '#pmolends';
                                    returnAnswer(answer,parentID,units);
                                };
                            }
            break;
        case "ratio":
            calcFunction = 	function (e){
                                clearAnswerWarningField(e);
                                var inputFields = $('form').find('input');
                                if (validateInputValues(inputFields)) {
                                    var answer = parseFloat($("#ins").val()) * parseFloat($('#vectorng').val()) / parseFloat($('#vectorkb').val());
                                    answer = biomathRounding(answer);
                                    var units = "ng insert for 1:1 ratio";
                                    var parentID = '#ratio';
                                    returnAnswer(answer,parentID,units);
                                };
                            }
            break;
        case "odConvert":
            calcFunction = 	function (e){
                                
                                clearAnswerWarningField(e);
                                var inputFields = $('form').find('input');
                                if (validateInputValues(inputFields)) {
                                    if ($('#od').val() < 0.1 || $('#od').val() >1){
                                        $('#odConvertMessage').html("Most spectrophotometers are not very accurate below 0.1 OD units or above 1.0 OD units. You might want to consider measuring a different dilution of your sample.").show();
                                    }
                                    var answer = parseFloat($('#od').val()) * parseFloat($('#naType option:selected').val());
                                    answer = biomathRounding(answer);
                                    var units = "&micro;g/ml";
                                    var parentID = '#odConvert';
                                    returnAnswer(answer,parentID,units);                          
                                }
                            }
            break;
        case "kdapmolug":
            calcFunction =  function (e){
                                clearAnswerWarningField(e);

                                var inputFields = $('form').find('input').not(':hidden,input[type=hidden]');
                                if (validateInputValues(inputFields)) {
                                    var kda=0;
                                    var pmol=0;
                                    var ug=0;
                                    var selected = $('#kdapmolugSelect option:selected').val();
                                    
                        
                                    if (selected == '1'){
                                        kda= parseFloat($('#kdaProtein').val());
                                        pmol= parseFloat($('#pmolProtein').val());  
                                        ug= Math.round(kda * pmol) * 0.001;
                                        var answer = biomathRounding(ug);
                                        var units = "&micro;g of protein";
                                        var parentID = '#kdapmolug';
                                        returnAnswer(answer, parentID, units);
                                    }
                                    if (selected == '2'){
                                        kda= parseFloat($('#kdaProtein').val());
                                        ug = parseFloat($('#ugProtein').val());     
                                        pmol = (1000 * ug) /kda;
                                        var answer = biomathRounding(pmol);
                                        var units = "pmol of protein";
                                        var parentID = "#kdapmolug";
                                        returnAnswer(answer, parentID, units);
                                        $('#pmolProtein').closest('.control-group').removeClass('error');
                                    }
                                    if (selected == '0'){
                                        pmol= parseFloat($('#pmolProtein').val());
                                        ug = parseFloat($('#ugProtein').val());     
                                        kda = (1000 * ug) /pmol;
                                        var answer = biomathRounding(kda);
                                        var units = "kDa";
                                        var parentID = "#kdapmolug";
                                        returnAnswer(answer, parentID, units);
                                        $('#kdaProtein').closest('.control-group').removeClass('error');
                                    }
                                }
                                
                            }
            //Custom clear function
            clearFunction = function (e){
                                clearAnswerWarningField(e);
                                $("#kdaProtein").val("");
                                $("#kdaProtein").closest('.control-group').removeClass('error');
                                $("#ugProtein").val("");
                                $("#ugProtein").closest('.control-group').removeClass('error');
                                $("#pmolProtein").val("");
                                $("#pmolProtein").closest('.control-group').removeClass('error');
                            }
            break;
        case "dnaprotein":
            calcFunction =  function (e){
                                clearAnswerWarningField(e);

                                var inputFields = [mostRecentDnaproteinField];

                                if (validateInputValues(inputFields)) {
                                    var kda = 0;
                                    var aminoAcid = 0;
                                    var dnaSize = 0;
                                    
                                    if ('kda' === mostRecentDnaproteinField.attr('id')){
                                        kda = parseFloat($('#kda').val());
                                        aminoAcid = Math.round(kda / 0.11);
                                        dnaSize = aminoAcid * 3;
                                        $('#aminoAcid').val(biomathRounding(aminoAcid));
                                        $('#dnaSize').val(biomathRounding(dnaSize));
                                        return;
                                    }
                                    if ('aminoAcid' === mostRecentDnaproteinField.attr('id')){
                                        aminoAcid = parseFloat($('#aminoAcid').val());
                                        kda = aminoAcid * 0.11;
                                        dnaSize = aminoAcid * 3;
                                        $('#kda').val(biomathRounding(kda));
                                        $('#dnaSize').val(biomathRounding(dnaSize));
                                        return;
                                    }
                                    if ('dnaSize' === mostRecentDnaproteinField.attr('id')) {
                                        dnaSize = parseFloat($('#dnaSize').val());
                                        aminoAcid = Math.floor(dnaSize / 3);
                                        kda = aminoAcid * 0.11;
                                        $('#aminoAcid').val(biomathRounding(aminoAcid));
                                        $('#kda').val(biomathRounding(kda));
                                        return;
                                    }
                                };

                                
                            }
            //Custom clear function
            clearFunction = function (e){
                                clearAnswerWarningField(e);
                                $('#dnaproteinForm input').val('');
                            }
            break;
        case "dilution":
            calcFunction =  function (e){
                                clearAnswerWarningField(e);
                                var inputFields = $('form').find('input');
                                if (validateInputValues(inputFields)) {
                                    //Units for concentration of solution used by user  Actually records the power of 10 from standard unit
                                    var intStockPower = $('[name="startunit"] option:selected').val();
                                    //Units for concentration of solution desired by user  Actually records the power of 10 from standard unit          
                                    var intFinalPower = $('[name="finalunit"] option:selected').val();  
                                    //Units for volume of solution desired by user  Actually records the power of 10 from standard unit     
                                    var intVolumePower = $('[name="volunit"] option:selected').val();
                                    //Absolute concentration of stock solution      
                                    var sngStockConcentration = parseFloat($('#startconc').val()) * Math.pow(10,intStockPower);
                                    //Absolute concentration of final solution  
                                    var sngFinalConcentration = parseFloat($('#finalconc').val()) * Math.pow(10,intFinalPower);
                                    //Absolute volume of final solution 
                                    var volume = parseFloat($('#finalvol').val()) * Math.pow(10,intVolumePower);
                                    //Dilution factor needed to make their final solution                   
                                    var sngDilution = sngStockConcentration/sngFinalConcentration;
                                    //Absolute volume of stock needed               
                                    var sngNeededVolume = volume/sngDilution;   
                                    //'Reported volume of stock needed. Will have associated units.     
                                    var sngResultVolume = 0;
                                    //Name of units to be reported          
                                    var strResultUnit = "";
                                    //Extra warning string
                                    var strWarning = "";
                                    
                                    function finalDilutionAnswer(strWarning,sngResultVolume,strResultUnit){
                                        $('#dilutionAnswer').html("Add " + sngResultVolume + strResultUnit + " of stock solution.").addClass('answer').show();
                                        if (strWarning !=='') {
                                            alert(strWarning);
                                        };
                                    }

                                    if (checkForNaNorInfinity(sngNeededVolume)) {
                                        $('#dilutionAnswer').html("Please check for an invalid entry.").removeClass('answer').show(); 
                                        return;
                                    };
                                
                                    if (sngFinalConcentration > sngStockConcentration){
                                        alert("Your stock solution is less concentrated than your final solution. You would need to concentrate your stock solution " + 1/sngDilution + " fold.");
                                    }
                                    else{ 
                                        if (sngNeededVolume >= 1){
                                            sngResultVolume = sngNeededVolume;
                                            sngResultVolume = Math.round(sngResultVolume*1000)/1000;
                                            strResultUnit=" L"; 
                                            finalDilutionAnswer(strWarning,sngResultVolume,strResultUnit);
                                        }
                                        else {
                                            if (sngNeededVolume >=Math.pow(10,-3)){
                                                sngResultVolume=sngNeededVolume * Math.pow(10,3);
                                                sngResultVolume=Math.round(sngResultVolume*100)/100;
                                                strResultUnit="ml";
                                                finalDilutionAnswer(strWarning,sngResultVolume,strResultUnit);
                                            }
                                            else{ 
                                                if (sngNeededVolume >= Math.pow(10,-6)){
                                                    strResultUnit="&micro;l";
                                                    sngResultVolume=sngNeededVolume * Math.pow(10,6);
                                                    sngResultVolume=Math.round(sngResultVolume*10)/10;
                                                    finalDilutionAnswer(strWarning,sngResultVolume,strResultUnit);
                                                }
                                                    
                                                else{ 
                                                    if (sngNeededVolume >=Math.pow(10,-9)){
                                                        sngResultVolume=sngNeededVolume * Math.pow(10,9);
                                                        sngResultVolume=Math.round(sngResultVolume);
                                                        strResultUnit="nl";
                                                        strWarning="It is not generally possible to pipette accurately at these volumes. You may want to consider diluting your stock solution. ";
                                                        finalDilutionAnswer(strWarning,sngResultVolume,strResultUnit);
                                                    }
                                                    else {
                                                        if (sngNeededVolume >=Math.pow(10,-12)){                            
                                                            sngResultVolume=sngNeededVolume * Math.pow(10,12); 
                                                            sngResultVolume=Math.round(sngResultVolume);
                                                            strResultUnit="pl";
                                                            strWarning="It is not generally possible to pipette accurately at these volumes. You may want to consider diluting your stock solution. ";
                                                            finalDilutionAnswer(strWarning,sngResultVolume,strResultUnit);
                                                        }
                                                        else {
                                                            
                                                            alert("The volume needed is less than one picoliter. Please make sure that you have entered the correct information. You will need to make a significant dilution of your stock solution to accurately make your desired final solution.");
                                                        }
                                                    }
                                                }
                                            }           
                                        }
                                    }
                                }
                            }
            break;
        case "molarity":
            calcFunction = function (e){
                                clearAnswerWarningField(e);
                                var inputFields = $('form').find('input');
                                if (validateInputValues(inputFields)) {
                                    var  sngVolume = parseFloat($('#volume').val());    //Volume of solution desired by user
                                    var  intVolume_Units = $('[name="volumeUnits"] option:selected').val(); //Units for volume of solution desired by user  Actually records the power of 10 from standard unit
                                    var  sngMolWt = parseFloat($('#mw').val());         //Molecular weight of compound entered by user
                                    var  sngMolarity = parseFloat($('#concFinal').val());       //Molarity of of solution desired by user
                                    var  intMolarity_Units = $('[name="molarityUnits"] option:selected').val(); //Units for molarity of solution desired by user. Actually records the power of 10 from standard unit
                                    var  sngGrams_Needed;   //Absolute number of grams needed for the results
                                    var  sngResultGrams;        //Number of grams to report. Will also have units associated
                                    var  strResultUnits;        //Units to report
                                //  var  strVolUnitName;        //text for volume units
                                //  var  strMolUnitName;        //text for concentration units
                                //  var strWarning = "";        //additional warning
                                    sngGrams_Needed=(sngVolume * Math.pow(10,intVolume_Units) * sngMolarity * Math.pow(10,intMolarity_Units) * sngMolWt);
                                    
                                    if (checkForNaNorInfinity(sngGrams_Needed)) {
                                        $('#molarityAnswer').html("Please check for an invalid entry.").removeClass('answer').show(); 
                                    };

                                    strResultUnits=" grams";
                                    if (sngGrams_Needed === 1) {
                                        strResultUnits = " gram";
                                    }

                                    if  (sngGrams_Needed >= 1){
                                        sngResultGrams=sngGrams_Needed;
                                        sngResultGrams=Math.round(sngResultGrams*1000)/1000;

                                        $('#molarityAnswer').html("You will need " + sngResultGrams + strResultUnits + ".").addClass('answer').show();   
                                    }
                                    else if(sngGrams_Needed >= Math.pow(10,-3) && sngGrams_Needed <1){          
                                        sngResultGrams=sngGrams_Needed * Math.pow(10,3);
                                        sngResultGrams=Math.round(sngResultGrams*100)/100;
                                        strResultUnits=" milligrams";
                                        $('#molarityAnswer').html("You will need " + sngResultGrams + strResultUnits + ".").addClass('answer').show();       
                                    }
                                    else if(sngGrams_Needed >= Math.pow(10,-6) && sngGrams_Needed <Math.pow(10,-3)){
                                        sngResultGrams=sngGrams_Needed * Math.pow(10,6);
                                        sngResultGrams=Math.round(sngResultGrams*10)/10;
                                        strResultUnits=" micrograms";
                                        alert("It is not generally possible to accurately weight this amount. You may want to consider making a more concentrated stock solution and then diluting to your desired final concentration.");
                                        $('#molarityAnswer').html("You will need " + sngResultGrams + strResultUnits + ".").addClass('answer').show();   
                                    }
                                    else if(sngGrams_Needed < Math.pow(10,-6)){
                                        alert("The amount needed is less than one microgram. Please make sure that you have entered the correct information. You will need to make a significantly more concentrated stock solution and then dilute to make your desired final solution.");
                                    }
                                }
                            }
            break;
        case "temp":
            calcFunction =  function (e){


                                clearAnswerWarningField(e);

                                var inputFields = [mostRecentTempField];

                                if (validateInputValues(inputFields)) {
                                
                                    var f = parseFloat($('#f').val());
                                    var c = parseFloat($('#c').val());
                                    var k = parseFloat($('#k').val());
                                    
                                    if ('f' === mostRecentTempField.attr('id')){
                                        c= 5/9 *(f-32);
                                        k= parseFloat(c) + 273.16 ;
                                    }
                                    if ('c' === mostRecentTempField.attr('id')){
                                        f= 9/5 *c +32;
                                        k= parseFloat(c) + 273.16;
                                    }
                                    if ('k' === mostRecentTempField.attr('id')){
                                        c= k -  273.16;
                                        f= 9/5 *c +32;
                                    }
                                    if(k < 0){
                                        $('#kWarning').html("This value is below absolute zero.");
                                        $('#kWarning').closest('.calcField').removeClass('warning info').addClass('error');
                                    }else{
                                        $('#kWarning').closest('.calcField').removeClass('error');
                                    $('#f').val((Math.round(f*100))/100);
                                    $('#k').val((Math.round(k*100))/100);
                                    $('#c').val((Math.round(c*100))/100);
                                    }
                                };
                                
                                
                            }	
            //Custom clear function
            clearFunction = function (e){
                                clearAnswerWarningField(e);
                                $('#tempForm input').val('');
                            }
            break;
        case "tm":
            calcFunction =  function (e){
                                clearAnswerWarningField(e);
                                function countbase(oligo, base){//counts the number of times a specific base appears in the sequence
                                    var number= 0;
                                    var other = 0;
                                    for (var i = 0;  i < oligo.length;  i++)
                                    {
                                        var ch = oligo.charAt(i);
                                        if (ch == base)
                                        {
                                            number++;
                                        }
                                    }
                                    return number;
                                }
                                   
                                function gcContent(realseq){
                                    var G=countbase(realseq, "G");
                                    var A=countbase(realseq, "A");
                                    var T=countbase(realseq, "T");
                                    var C=countbase(realseq, "C");
                                    var GC=100*((G+C)/(G+A+T+C));
                                    var gcAnswer=Math.round((100*GC)/100);
                                    return gcAnswer;
                                }
                                
                                function cleanseq(seq){ //changes all bases to uppercase for easier manipulation and removes any nonbase characters.
                                    
                                    var newseq = "";
                                    seq=seq.toUpperCase();
                                    var counter=1;
                                    for (var i = 0; i < seq.length; i++)
                                    {	
                                        if ((seq.charAt(i)=="G")||(seq.charAt(i)== "A")||(seq.charAt(i)== "T")||(seq.charAt(i)== "C"))
                                        {
                                            newseq=newseq+seq.charAt(i);
                                        }else{
                                            counter++;
                                        }
                                    }
                                    if(counter>1){	
                                        alert("Characters other than G, C, T and A will be removed from your sequence.");
                                    }
                                    return newseq;
                                }
                                
                                // base stacking calculations. 
                                function stack(seq,salt,primer,mg){
                                   var tmAnswer;
                                   var dh=0;
                                   var ds=0;
                                   var K=(primer/2) * Math.pow(10,-9); //converts from nanomolar to Molar. Note this ignores the contribution of the target since this is << than primer concentration.
                                   
                               //Salt corrections calculations. Used to be saltcorrections().
                               //changes to ds dependant on the salt concentration & sequence length
                                   
                                   salt+=(mg/Math.pow(10,3) * 140);//convert to moles and then adjust for greater stabilizing effects of Mg compared to Na or K. See von Ahsen et al 1999	
                                    var deltas=0;
                                   deltas+=0.368 * (seq.length-1)* Math.log(salt);//This comes from von Ahsen et al 1999
                                   ds+=deltas;
                               //Used to be terminalcorrections()
                               //helix initiation corrections from Santalucia 1998 & Allawi & Santalucia 1997
                                   var deltah=0;
                                   var deltas=0;
                                   if ((seq.charAt(0)=="G")||(seq.charAt(0)=="C"))
                                   {
                                       deltah+=0.1;
                                       deltas+=-2.8;
                                   }
                                   if((seq.charAt(0)=="A")||(seq.charAt(0)=="T"))
                                   {
                                       deltah+=2.3;
                                       deltas+=4.1;
                                   }
                                   
                                   if ((seq.charAt(seq.length-1)=="G")||(seq.charAt(seq.length-1)=="C"))
                                   {
                                       deltah+=0.1;
                                       deltas+=-2.8;
                                   }
                                   if((seq.charAt(seq.length-1)=="A")||(seq.charAt(seq.length-1)=="T"))
                                   {
                                       deltah+=2.3;
                                       deltas+=4.1;
                                   }
                                   dh+=deltah;
                                   ds+=deltas;
                               //end terminalcorrections()
                                   
                                   for (var i = 0; i < seq.length; i++)//adds up dh and ds for each 2 base combination. dh is in kcal/mol. ds is in cal/Kelvin/mol
                                   {	
                                       if (seq.charAt(i)=="G")
                                       {
                                           if (seq.charAt(i+1)=="G")
                                           {
                                               dh+=-8;
                                               ds+=-19.9;
                                           }
                                           if (seq.charAt(i+1)=="A")
                                           {
                                               dh+=-8.2;
                                               ds+=-22.2;
                                           }
                                           if (seq.charAt(i+1)=="T")
                                           {
                                               dh+=-8.4;
                                               ds+=-22.4;
                                           }
                                           if (seq.charAt(i+1)=="C")
                                               //These values where fixed on 4/23/2008. They were dh = -10.6 & ds = -27.2
                                               //The new values have been double check against Santalucia 1998 
                                           {
                                               dh+=-9.8;
                                               ds+=-24.4;
                                           }
                                       }
                                       if (seq.charAt(i)=="A")
                                       {
                                           if (seq.charAt(i+1)=="G")
                                           {
                                               dh+=-7.8;
                                               ds+=-21;
                                           }
                                           if (seq.charAt(i+1)=="A")
                                           {
                                               dh+=-7.9;
                                               ds+=-22.2;
                                           }
                                           if (seq.charAt(i+1)=="T")
                                           {
                                               dh+=-7.2;
                                               ds+=-20.4;
                                           }
                                           if (seq.charAt(i+1)=="C")
                                           {
                                               dh+=-8.4;
                                               ds+=-22.4;
                                           }
                                       }
                                       if (seq.charAt(i)=="T")
                                       {
                                           if (seq.charAt(i+1)=="G")
                                           {
                                               dh+=-8.5;
                                               ds+=-22.7;
                                           }
                                           if (seq.charAt(i+1)=="A")
                                           {
                                               dh+=-7.2;
                                               ds+=-21.3;
                                           }
                                           if (seq.charAt(i+1)=="T")
                                           {
                                               dh+=-7.9;
                                               ds+=-22.2;
                                           }
                                           if (seq.charAt(i+1)=="C")
                                           {
                                               dh+=-8.2;
                                               ds+=-22.2;
                                           }
                                       }
                                       if (seq.charAt(i)=="C")
                                       {
                                           if (seq.charAt(i+1)=="G")
                                           {
                                               dh+=-10.6;
                                               ds+=-27.2;
                                           }
                                           if (seq.charAt(i+1)=="A")
                                           {
                                               dh+=-8.5;
                                               ds+=-22.7;
                                           }
                                           if (seq.charAt(i+1)=="T")
                                           {
                                               dh+=-7.8;
                                               ds+=-21;
                                           }
                                           if (seq.charAt(i+1)=="C")
                                           {
                                               dh+=-8;
                                               ds+=-19.9;
                                           }
                                       }
                                   }
                                   tmAnswer = Math.round(((1000* dh)/(ds+(1.987 * Math.log(K))))-273.15);
                                   return tmAnswer;				
                               }
                                
                                var mg=1.5;
                                var seq = $('#sequence').val();
                                var realseq = cleanseq(seq);//remove non-bases from the sequence
                                $('#sequence').val(realseq); //show the sequence that will actually be used for calculations
                                 var salt=50/Math.pow(10,3);//convert to molar Math.pow(10,-3))
                                 var primer=200;
                                 //This if statement prevents calculator from returning NaN if a required field was empty
                                // if(realseq&&salt&&primer){
                                var tmAnswer=stack(realseq,salt,primer,mg);
                                var gcAnswer=gcContent(realseq);

                                if(checkForNaNorInfinity(gcAnswer)){
                                    $('#tmAnswer').html("Please check for an invalid entry.").removeClass('answer').show();
                                }else{
                                    $('#tmAnswer').html("%GC: " + gcAnswer + " and Tm: " + tmAnswer +"&deg;C").addClass('answer').show();
                                }
                            }
                        clearFunction = function (e){
                                clearAnswerWarningField(e);
                                $('#tmForm input').val('');
                                $('#tmForm select').prop('selectedIndex', 0);
                            }
            break;
        default:
            break;
    }
        //Create the calculator object & add to calculators array
        var calculator = new Calculator(calculatorIdArray[i], calculatorNameArray[i], number, calcFunction, clearFunction);
        calculators[i] = calculator;
    }
}



function Calculator(id, name, number, calcFunction, clearFunction)
{
    //Calculator Name/html ID tag
    this.id = id;
    this.name = name;
    //Calculator sequential number used for formula image filename
    this.number = number;
    //Calculation function
    this.calcFunction = calcFunction;
    //Clear function
    this.clearFunction = clearFunction;
    
    this.bindListHandler();
}

Calculator.prototype.bindListHandler = function()
{
    var calc = this;
    
    $('#' + calc.id + "List").bind('click', function() {
        calc.loadCalculator();
    });
}

Calculator.prototype.loadCalculator = function()
{
    var calc = this;

    $('#centerPanel').load("html_imports/" + calc.id + ".html", function() {
        if (calc.id === 'kdapmolug') {
            kdapmolugSelectSetup();
        };
        $('#calculateButton').bind("click", calc.calcFunction).bind('touchstart', function(){$('#calculateButton').addClass('buttonTapped'); }).bind('touchend', function(){$('#calculateButton').removeClass('buttonTapped');});
        $('#clearButton').bind("click", calc.clearFunction).bind('touchstart', function(){$('#clearButton').addClass('buttonTapped'); }).bind('touchend', function(){$('#clearButton').removeClass('buttonTapped');});

        loadFieldDisabler(calc.id);
    });
    $('#formulaDiv').load("html_imports/" + calc.id + "_calc.html", function() {
        var formulaDiv = document.getElementById('formulaDiv');
        MathJax.Hub.Queue(["Typeset", MathJax.Hub, formulaDiv]);
    });
    $('#title').html(calc.name);

    if (!isInitialRequest) {
        _gaq.push(['_trackPageview', "/a/apps/biomath/index.html?calc=" + this.id]);
    } else {
        isInitialRequest = false;
    }
    


    loadCalculateClearButtons();
    
}

function kdapmolugSelectSetup() {
    $('#proteinSize').hide();
    $('#kdapmolugSelect').change(function() {
        var selected = $(this).find(":selected").val();
        
        switch (selected) {
            case '0':
                $('#proteinSize').closest('.control-group').hide();
                $('#ugProtein').closest('.control-group').show();
                $('#pmolProtein').closest('.control-group').show();
                break;
            case '1':
                $('#proteinSize').closest('.control-group').show();
                $('#ugProtein').closest('.control-group').hide();
                $('#pmolProtein').closest('.control-group').show();
                break;
            case '2':
                $('#proteinSize').closest('.control-group').show();
                $('#ugProtein').closest('.control-group').show();
                $('#pmolProtein').closest('.control-group').hide();
                break;
            default:
                break;
        }
        $('#clearButton').click();
    });
}

function loadCalculateClearButtons()
{
    //Resets Calculate and Clear buttons when new values are being entered
	$('input:not(.noClearOnFocus)').focusin(resetButton); 
	$('select').change(resetButton);
}

function loadFieldDisabler(id)
{

    if (id == "dnaprotein") {
        console.log('dnaprotein');

        $('#dnaproteinForm input').each(function(){
            console.log('hello i am a number input');
        });
        
        $('#dnaproteinForm input').bind('keyup', function(){

            console.log('keyup');
            
            if ($(this) === mostRecentDnaproteinField) {
                return;
            } else {
                mostRecentDnaproteinField  = $(this);
            }

            if ($(this).attr('id') === 'dnaSize') {
                $('#aminoAcid').val('');
                $('#kda').val('');
            };
            if ($(this).attr('id') === 'aminoAcid') {
                $('#dnaSize').val('');
                $('#kda').val('');
            };
            if ($(this).attr('id') === 'kda') {
                $('#dnaSize').val('');
                $('#aminoAcid').val('');
            };
		});
    } else if (id == "temp"){
        console.log('temp');

        $('#tempForm input').bind('keyup', function(){
            console.log('keyup');

            if ($(this) === mostRecentTempField) {
                return;
            } else {
                mostRecentTempField = $(this);
            }

            if ($(this).attr('id') === 'f') {
                $('#c').val('');
                $('#k').val('');
            };
            if ($(this).attr('id') === 'c') {
                $('#f').val('');
                $('#k').val('');
            };
            if ($(this).attr('id') === 'k') {
                $('#c').val('');
                $('#f').val('');
            };
        });
    } else if (id == "tm") {
        
        $('#primers').bind('change', function(){
			document.getElementById("sequence").value=this.value;
			document.getElementById("primers").focus();	
		});
    }
}
 
/** 
 * UTILITY FUNCTIONS - Used by multiple calculators
 * */

//Returns the Calculate and Clear buttons
function resetButton() {
    var parentID = "#"+$(this).closest(".calculator").attr('id');
    $(parentID+"Answer").hide().html("").removeClass('answer');
    $(parentID+"Btn").show();
    $('.messageHolder').hide();			
}

//Calc answers are passed to this function. Warning is displayed or answer is returned to user.		
function returnAnswer(calcAnswer,parentID,units) { 
    if (checkForNaNorInfinity(calcAnswer)) {
        $(parentID+'Answer').html("Please check for an invalid entry.").show();
    } else {
        $(parentID+'Answer').html(insertCommas(calcAnswer)+units).addClass('answer').show();
    }
}

function insertCommas(rawAnswer) {
    if (rawAnswer.toString().indexOf('.') != -1) {
        var numberParts = rawAnswer.toString().split('.');
        var leftOfDecimal = numberParts[0].replace(/\B(?=(\d{3})+(?!\d))/g, ",");
        return (leftOfDecimal + "." + numberParts[1]); 
    } else {
        return rawAnswer.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
    }
}

function biomathRounding(answer){
    if(checkForNaNorInfinity(answer)){
        //display error
        return answer;
    }
    if(answer >= .1){
        return answer = Math.round(answer * 100) / 100;
    }
    if(answer >= .01 && answer < .1){
        return answer = Math.round(answer  *1000) / 1000;
    }
    
    if(answer >= .001 && answer < .01){
        return answer = Math.round(answer * 10000) / 10000;
    }
    if(answer >= .0001 && answer < .001){
        return answer = Math.round(answer * 100000) / 100000;   
    }
    if(answer < .0001){
        return answer = 0;  
    }
    return answer;
}

/**
 * Validation Code
 * */

function validateInputValues(inputFields) {

    var returnVal = true;


    $(inputFields).each(function(index){

        var fieldMinValue = $(this).attr('min');
        var warningFieldID = '#' + $(this).attr('id') + 'Warning';

        if (parseInt(fieldMinValue) === 1) {
            if (!positiveInteger($(this).attr('id'))) {
                $(warningFieldID).html("Please enter a positive integer.");
                $(warningFieldID).closest('.calcField').removeClass('warning info').addClass('error').show();
                returnVal = false;
            }
        };

        if (parseInt(fieldMinValue) === 0) {
            if (!positiveDecimal($(this).attr('id'))) {
                $(warningFieldID).html("Please enter a positive decimal.");
                $(warningFieldID).closest('.calcField').removeClass('warning info').addClass('error').show();
                returnVal = false;
            };
        };

        if (isNaN(fieldMinValue)) {
            if (decimal($(this).attr('id')) === false) {
                $(warningFieldID).html("Please check for an invalid entry.");
                $(warningFieldID).closest('.calcField').removeClass('warning info').addClass('error').show();
                returnVal = false;
            };
        };

        if ($(this).val() === '') {
            //empty field
            $(warningFieldID).html("Please check for a missing or invalid entry.");
            $(warningFieldID).closest('.calcField').removeClass('warning info').addClass('error').show();
            returnVal = false;
        }
    });
    return returnVal;
}

function decimal(fieldID){
    var validateThis = $('#' + fieldID).val();
        if(isNaN(validateThis)){
            return false;
        } else {
            return true;
        }
}
function positiveDecimal(fieldID){
    var validateThis = $('#' + fieldID).val();
    if(isNaN(validateThis) || parseFloat(validateThis) < 0){
        return false;
    } else {
        return true;
    }
}
function positiveInteger(fieldID){
    var validateThis = $('#' + fieldID).val();
    if(isNaN(validateThis) || parseFloat(validateThis) < 0 || parseFloat(validateThis) % 1 != 0){
        return false;
    } else {
        return true;
    }
} 

function checkForNaNorInfinity(result) {
    if(isNaN(result)||result===Infinity) {
        return true;
    } else {
        return false;
    }
}


