Application.Dna.service('PCR', ['Clotho', 'DNA', 'Digest', function(Clotho, DNA, Digest) {

    /**************
     Primer Design
     **************/

    /**
     * @description Get clamp primer using beginning of sequence
     * @param sequence
     * @param minTemp Temperature (°C) primer should be. Default 50°C
     * @returns {string} primer reverse complimentary to sequence, minimum 10nt, maximum 50nt
     */
    var generatePrimer = function generatePrimer(sequence, minTemp) {
        minTemp = !!minTemp ? minTemp : 50;

        var minLength = 10,
            maxLength = 50,
            startTest = Math.floor(minTemp / 3.5), //assumes GC_content < 75%
            index = (startTest > minLength ? startTest : minLength),
            primer = sequence.substring(0, index);

        while( DNA.melting_temp_basic(primer) < minTemp || index > maxLength){
            index = index +1;
            primer = sequence.substring(0, index);
        }

        return DNA.revcomp(primer);
    };

    /**************
     Annealing
     **************/

    /**
     * @description Finds indices where primer anneals to sequence, as given by Digest.findIndices() -- i.e. exact match without tail
     * @param sequence
     * @param primer
     * @returns {Object} In form {forward: Array, reverse: Array}} where each array is indices the primer anneals, empty if no matches.
     */
    var findAnnealFull = function findAnnealFull(sequence, primer) {

        //future - move to fuzzy search, account for tail on primer

        var forward = Digest.findIndices(sequence, primer, false);
        var reverse = Digest.findIndices(sequence, DNA.revcomp(primer), false);

        return {forward: forward, reverse: reverse}
    };


    var findAnnealFuzzy = function(sequence, primer) {
        //todo
    };

    var findAnnealSimple = function(template, primer) {
        //start from 3', go back until have unique match
        var searchFrag, searchReg, matches = {}, result;

        for(var initialBack = 8,
                start = primer.length - initialBack;
            searchFrag = primer.substring(start),
                searchReg = Digest.createRegex(searchFrag);
            start++)
        {

            //check forward
            matches.forward = template.match(searchReg);
            //check revcomp
            matches.reverse = DNA.revcomp(template).match(searchReg);

            if (!(matches.forward.length || matches.reverse.length)) {
                console.log('no *exact* matches found for length ' + start + ' from 3 prime end');
                result = false;
                break;
            } else if ((matches.forward.length + matches.reverse.length) == 1) {
                //return index of single match
                //todo - store forward or reverse?
                result = matches.forward.length ? template.search(searchReg) : DNA.revcomp(template).search(searchReg);
                break;
            }
        }

        return result;
    };

    var anneal = function (template, primer) {
        var method = findAnnealSimple;

        //todo

    };

    //e.g. PCA, given array of oligos, find which match which
    var annealArray = function (oligos) {

        var track = {};
        for (var i = 0; i < oligos.length; i++) {
            track[i] = {oligo : oligos[i]};

            //todo - find matches in oligos
        }

    };


    //TODO - write algorithm that runs the PCR, not as predictive

    /**************
     PCR algorithms
     ****************/


    /** verification **/
    /*
     multiple matches:
     - will work, multiple (unexpected?) products
     one match each: {p1_fwd,  p1_rev, p1_absent} x {p2_fwd,  p2_rev, p2_absent}
     - 3 cases (fails): 1 or both absent
     - 2 cases (fails): same direction
     - 2 cases: point toward each other
     p1 = fwd, p2 = rev (fwd < rev)  or p1 = rev, p2 = fwd (fwd < rev)
     - 2 cases: point away from each other
     same as above but (fwd > rev)
     special cases: ?
     */
    /**
     * @description Determines if primers will anneal and only once (Simple PCR)
     * @param sequence
     * @param primer1
     * @param primer2
     * @returns {string|boolean} true if no error, otherwise string with error
     */
    var verifyPrimers = function verifyPrimers(sequence, primer1, primer2) {

        var p1 = findAnnealFull(sequence, primer1),
            p2 = findAnnealFull(sequence, primer2);

        /** check zero matches **/

        if (p1.forward.length + p1.reverse.length < 1) {
            return "primer1 : no matches"
        }
        if (p2.forward.length + p2.reverse.length < 1) {
            return "primer2 : no matches"
        }

        /** check multiple matches **/

        if (p1.forward.length + p1.reverse.length > 1) {
            return "primer1 : multiple matches"
        }
        if (p2.forward.length + p2.reverse.length > 1) {
            return "primer2 : multiple matches"
        }

        /** check directions */

        if ((p1.forward.length && p2.forward.length) || (p1.reverse.length && p2.reverse.length))
            return "primers point same direction";

        return true;

    };



    // see https://www.ncbi.nlm.nih.gov/tools/epcr/

    //wrapper function, currently only handles PCR and EIPCR
    var predict = function predict(sequence, primers) {

        //console.log(sequence, primers);

        if (primers.length != 2)
            return "Can only handle having two primers right now";

        var primer1 = primers[0],
            primer2 = primers[1];

        var verify = verifyPrimers(sequence, primer1, primer2);
        if (verify !== true)
            return verify;


        //future - not DRY
        var p1 = findAnnealFull(sequence, primer1),
            p2 = findAnnealFull(sequence, primer2);


        /** orient primers **/
        //todo - break out

        var p1F = !!p1.forward.length;
        var p2F = !!p2.forward.length;
        var p1pos = p1F ? +p1.forward[0] : +p1.reverse[0];
        var p2pos = p2F ? +p2.forward[0] : +p2.reverse[0];



        //pass to protocol


        //p1 forward, p2 reverse
        if (p1F) {
            //normal
            if (p1pos < p2pos) {
                return PCR(sequence, p1pos, (p2pos + primer2.length));
            }
            //eipcr
            else {
                return EIPCR(sequence, p1pos, (p2pos+ primer2.length));
            }
        }
        //p2 forward, p1 reverse
        else {
            //normal
            if (p1pos > p2pos) {
                return PCR(sequence, p2pos, (p1pos + primer1.length));
            }
            //eipcr
            else {
                return EIPCR(sequence, p2pos, (p1pos + primer1.length));
            }

        }
    };



    /** extension **/
    //note- currently logic for this functions exists in wrapper function

    var PCR = function PCR(sequence, forwardPrimerPos, reversePrimerPos) {
        //console.log(arguments);
        return sequence.substring(forwardPrimerPos, reversePrimerPos);
    };

    //aka Overlap Extension
    //note - currently, very simple, only two overlapping primers with defined overlap length
    var SOEing = function () {

    };

    var EIPCR = function EIPCR(sequence, forwardPrimerPos, reversePrimerPos) {
        //console.log(arguments);
        return sequence.substring(forwardPrimerPos) + sequence.substring(0, reversePrimerPos);
    };

    //30-150 bp, 2 primers w/ homology region (~20bp) and extend
    var wobble = function wobble(primer1, primer2, overlapLength) {
        //future - calculate overlapLength if not given

        return primer1 + (DNA.revcomp(primer2)).substring(overlapLength);
    };

    var PCA = function PCA() {

    };

    var RTPCR = function RTPCR() {

    };

    var Gibson = function Gibson() {
        //adds ligase (remove nicks)
        //5' exo nuclease
    };

    var SLIC = function SLIC() {
        //3' exo nuclease

    };


    /*************
     Simple Enzyme Functions
     *************/

    /**
     * @description
     * @param sequence
     *
     * @example nnnnnn -> nnnnnn
     * @example nn_nnnn^ -> nn
     * @example nnnn_nnnnnnn -> nnnn
     * todo - handle beyond only sticky ends
     */
    var exonuclease35 = function(sequence) {

        if (sequence.indexOf('_') < 0)
            return sequence;

        var regex = /(.*?)(_.+?\^?.*)/gi,
            matches = regex.exec(sequence),
            removed = matches[2],
            remaining = matches[1];

        return remaining;
    };

    /**
     * @description
     * @param sequence
     *
     * @example
     * note - also handle beyond only sticky ends
     */
    var exonuclease53 = function(sequence) {

    };

    /**
     * @description
     * @param sequence
     *
     * @example nnnnnn -> nnnnnn
     * @example nn^nnnn -> ___________
     * @example nn^nnnn_ -> nnnnnn
     * @example nn_nnnn^nn -> ______________ exonuclease activity?
     *
     * todo - determine above, also handle beyond only sticky ends
     */
    var polymerase53 = function(sequence) {

    };




    /**************
     Alignment
     **************/


    /**************
     Prediction
     **************/




    /**************
     Construction - move to own module
     **************/


    var ligate = function ligate(fragments) {

        var overhangs = [],
            blunts = [];

        for (var i = 0; i < fragments.length; i++) {
            var frag = fragments[i],
                ends = Digest.findOverhangs(frag);

            for (var j = 0; j < ends.length; j++) {
                if (ends[i].isBlunt) {
                    blunts.push({fragment : frag, ends : ends})
                } else {
                    overhangs.push({fragment : frag, ends : ends})
                }
            }

        }

        //todo
        //note - blunt ends can ligate in multiple directions

    };


    return {
        predict : predict

    }
}]);

'use strict';

Application.Dna.directive('pcrPredict', ['PCR', function(PCR) {

    return {
        restrict: 'A',
        require: 'ngModel',
        scope: {
            backbone: '=ngModel',
            primers: '='
        },
        link: function pcrPredictLink(scope, element, attr, ngModel) {
            ngModel.$render = function() {
                process();
            };

            element.bind('blur keyup change', function() {
                scope.$apply(process());
            });

            function process () {
                element.text(PCR.predict(scope.backbone, scope.primers));
            }
        }
    };
}]);