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
    var findAnnealAllExact = function findAnnealAllExact(sequence, primer) {

        //future - move to fuzzy search, account for tail on primer

        var forward = Digest.findIndices(sequence, primer, false);
        var reverse = Digest.findIndices(sequence, DNA.revcomp(primer), false);

        return {forward: forward, reverse: reverse}
    };


    var findAnnealFuzzy = function(sequence, primer) {
        //todo
    };

    //todo - return array of objects, each with index and if forward and has overhang
    var findAnnealSimple = function(template, primer) {
        //start from 3', go back until have unique match
        var searchFrag, searchReg, result = null;

        for(var initialBack = 8,
                start = primer.length - initialBack,
                matches = {};
            start > 0,
                searchFrag = primer.substring(start),
                searchReg = Digest.createRegex(searchFrag);
            --start)
        {
            console.log(start);

            //check forward
            matches.forward = template.match(searchReg) || [];
            //check revcomp
            matches.reverse = DNA.revcomp(template).match(searchReg) || [];

            console.log(matches, matches.forward.length, matches.reverse.length, matches.forward.length + matches.reverse.length);

            if (!matches.forward.length && !matches.reverse.length) {
                console.log('no *exact* matches found for length ' + start + ' from 3 prime end');
                //todo - step back one, check for single match
                result = false;
                break;
            } else if ((matches.forward.length + matches.reverse.length) == 1) {
                console.log('one match' );
                //todo - go back as far as possible
                //return index of single match
                result = matches.forward.length ? template.search(searchReg) : DNA.revcomp(template).search(searchReg);
                break;
            } else {
                console.log('else');
            }
        }

        return result;
    };

    //wrapper function...
    var anneal = function (template, primer, fuzzy) {
        fuzzy = !!fuzzy || true;

        //future - once write fuzzy, implement here
        return (fuzzy) ? findAnnealSimple(template, primer) : findAnnealSimple(template, primer);

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

        if (angular.isEmpty(primer1) || angular.isEmpty(primer2)) {
            return "a primer is not defined"
        }

        var p1 = findAnnealAllExact(sequence, primer1),
            p2 = findAnnealAllExact(sequence, primer2);

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
        var p1 = findAnnealAllExact(sequence, primer1),
            p2 = findAnnealAllExact(sequence, primer2);


        /** orient primers **/
        //todo - break out

        var p1pos = (!!p1.forward.length) ? +p1.forward[0] : +p1.reverse[0];
        var p2pos = (!!p2.forward.length) ? +p2.forward[0] : +p2.reverse[0];



        //pass to protocol

        //p1 forward, p2 reverse
        if (!!p1.forward.length) {
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
    //todo - fold into PCR function
    //todo - handle multiple lines
    var primerAlign = function(sequence, primers) {
        if (primers.length != 2)
            return "Can only handle having two primers right now";

        var primer1 = primers[0],
            primer2 = primers[1];

        var verify = verifyPrimers(sequence, primer1, primer2);
        if (verify !== true)
            return verify;


        //future - not DRY
        var p1 = findAnnealAllExact(sequence, primer1),
            p2 = findAnnealAllExact(sequence, primer2);


        /** orient primers **/
        //todo - break out

        var p1pos = (!!p1.forward.length) ? +p1.forward[0] : +p1.reverse[0];
        var p2pos = (!!p2.forward.length) ? +p2.forward[0] : +p2.reverse[0];

        //console.log(p1pos, p2pos, sequence, primers[0], primers[1]);

        //todo - don't assume p1 first

        var line1 = DNA.createRun(' ', p1pos);
            line1 += (!!p1.forward.length) ? primers[0] : DNA.revcomp(primers[0]);
            line1 += DNA.createRun(' ', (p2pos - p1pos - primers[0].length));
            line1 += (!!p2.forward.length) ? primers[1] : DNA.revcomp(primers[1]);
            line1 += DNA.createRun(' ', (sequence.length - line1.length));

        //note - handle line breaks outside function
        return (line1);
    };




    /**
     * @description Ligates two fragments. Optionally, shows alignment of ligation
     * @param fragments {Array} NOTE Currently only two
     * @param align {boolean} Whether to show alignment (both strands)
     * @param showHTML {boolean} Include HTML tags marking original fragments and complementarity region. Recommended when displaying as HTML. Default false.
     * @param showMarks {boolean} maintain cut marks. Not relevant if showHTML. Default false.
     * @returns {string}
     */
    var ligate = function(fragments, align, showHTML, showMarks) {
        var ends = parseFragmentEnds(fragments),
            blunts = ends.blunts,
            overhangs = ends.overhangs,
            fragPair;
        
        if (blunts.length > 1) {
            if (blunts.length == 2) {
                fragPair = orientFragmentsForJoin(blunts[0], blunts[1]);
            }
        }

        if (overhangs.length > 1) {
            if (overhangs.length == 2) {
                fragPair = orientFragmentsForJoin(overhangs[0], overhangs[1]);
            }
        }

        console.log(fragPair);

        if (!fragPair)
            return 'cut marks not defined';

        if (fragPair[0].overhang != DNA.revcomp(fragPair[1].overhang))
            return 'overhangs not complimentary (5\' orientation): ' + fragPair[0].overhang + ' /// ' + fragPair[1].overhang;


        var matchType = (fragPair[0].end.isBlunt) ? 'bluntmatch' : 'stickymatch';

        var finalText = '<ligate-frag>' + (fragPair[0].fragment).substring(0, fragPair[0].end.index) + '</ligate-frag>' +
            '<ligate-'+matchType+'>' + fragPair[0].end.match + '</ligate-'+matchType+'>' +
            (fragPair[1].fragment).substring(fragPair[1].end.index + fragPair[1].end.match.length);


        if (!!align) {
            var line2 = DNA.complement((fragPair[0].fragment).substring(0, fragPair[0].end.index)) +
                '<ligate-'+matchType+'>' + fragPair[1].end.match + '</ligate-'+matchType+'>' +
                '<ligate-frag>' + DNA.complement((fragPair[1].fragment).substring(fragPair[1].end.index + fragPair[1].end.match.length)) + '</ligate-frag>';

           finalText += "\n" + line2;
        }

        if (!showHTML) {
            finalText = finalText.replace(/(<([^>]+)>)/ig, '');
        }

        if (!showMarks && !showHTML) {
            finalText = Digest.removeMarks(finalText);
        }

        return finalText

    };


    /**************
     ligation
     **************/

    var parseFragmentEnds = function PCR_parseFragmentEnds (fragments) {
        var overhangs = [],
            blunts = [];

        for (var i = 0; i < fragments.length; i++) {
            var frag = fragments[i],
                ends = Digest.findOverhangs(frag);

            for (var j = 0; j < ends.length; j++) {
                //todo - handle non-terminal marks, maybe by cutting?
                if (!ends[j].terminal) {}

                if (ends[j].isBlunt) {
                    blunts.push({fragment : frag, end : ends[j], overhang : "|"})
                } else {
                    overhangs.push({fragment : frag, end : ends[j], overhang : ends[j][3]})
                }
            }
        }

        return {blunts: blunts, overhangs : overhangs};
    };


    var reverseFragment = function (frag) {
        frag.fragment = DNA.revcomp(frag.fragment);
        frag.end = Digest.findOverhangs(frag.fragment)[0]; // assumes only 1 end
        if (!frag.end.isBlunt) frag.overhang = frag.end[3];
        return frag;
    };

    // orients all fragments
    // fiveprime - default true: so showing 5' overhang (...^...._...)
    var orientFragments = function (fragments, fiveprime) {
        fiveprime = fiveprime || true;
        for (var i = 0; i < fragments.length; i++) {
            if (fragments[i].end.is3prime == fiveprime)
                fragments[i] = reverseFragment(fragments[i]);
        }
        return fragments;
    };

    //orient so frag1 comes first (i.e., frag1 end is at end) and have same direction overhang
    //todo - handle non-terminal ends
    var orientFragmentsForJoin = function (frag1, frag2) {
        //first process frag1
        if (frag1.end.index == 0 || frag1.end.index < (frag1.fragment.length/2)) {
            frag1 = reverseFragment(frag1);
        }
        //first, orient so longer section to right
        if (frag2.end.index > frag2.fragment.length/2) {
            frag2 = reverseFragment(frag2);
        }
        //???????
        if (frag2.overhang != DNA.revcomp(frag1.overhang)) {
            frag2 = reverseFragment(frag2);
        }

        return [frag1, frag2];
    };

    /**************
     Prediction
     **************/




    /**************
     Construction - own module
     **************/


    return {
        predict : predict,

        anneal : anneal,
        findAnnealSimple : findAnnealSimple,

        primerAlign : primerAlign,

        parseFragmentEnds : parseFragmentEnds,
        reverseFragment : reverseFragment,
        orientFragments : orientFragments,
        orientFragmentsForJoin : orientFragmentsForJoin,

        exonuclease35 : exonuclease35,
        exonuclease53 : exonuclease53,
        polymerase53 : polymerase53,

        ligate : ligate

    }
}]);

Application.Dna.directive('pcrPredict', ['PCR', 'Digest', 'DNA', function(PCR, Digest, DNA) {

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

            scope.$watch(function() {
                return scope.primers[0] + scope.primers[1]
            }, function() {
                process();
            });

            function process () {
                element.text(PCR.predict(scope.backbone, scope.primers));
            }
        }
    };
}]);

Application.Dna.directive('pcrAlign', ['PCR', 'Digest', 'DNA', '$filter', function(PCR, Digest, DNA, $filter) {
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

            scope.$watch('primers', process, true);

            function process () {
                var alignment = PCR.primerAlign(scope.backbone, scope.primers);

                //console.log(alignment);

                //todo -- pass in width of element to breaklines
                var charNum = 57;

                alignment = $filter('breakLines')(alignment, charNum, "*").split('*');

                var backboneText = $filter('breakLines')(scope.backbone, charNum, "*").split('*');

                //console.log(alignment, backboneText);

                var finalText = "";
                for (var i = 0; i < backboneText.length; i++) {
                    finalText += alignment[i] + "\n" + backboneText[i] + "\n";
                }

                //console.log(finalText);
                element.html(finalText)
            }



        }
    }
}]);

Application.Dna.directive('ligateAlign', ['PCR', 'Digest', 'DNA', '$compile', '$filter', function(PCR, Digest, DNA, $compile, $filter) {
    return {
        restrict: 'A',
        require: 'ngModel',
        scope: {
            fragments: '=ngModel'
        },
        link: function pcrPredictLink(scope, element, attr, ngModel) {
            ngModel.$render = function() {
                process();
            };

            scope.$watch('fragments', function() {
                process();
            }, true);

            function process () {
                var alignment = PCR.ligate(scope.fragments, true, true);
                console.log(alignment);

                alignment = $filter('DigestCuts')(alignment, true);

                element.html($compile('<span>' + alignment + '</span>')(scope))
            }

        }
    }
}]);

Application.Dna.directive('ligateFrag', function() {
    return {
        restrict: 'EA',
        replace: false,
        transclude:true,
        template: '<span tooltip="Initial Fragment" tooltip-placement="mouse" tooltip-animation="false" tooltip-append-to-body="true" ng-transclude></span>',
        link: function(scope, element, attrs) {
            element.css('color', '#faa');
        }
    }
});

Application.Dna.directive('ligateStickymatch', function() {
    return {
        restrict: 'EA',
        replace: false,
        transclude:true,
        template: '<span tooltip="Sticky-end Complementary Region" tooltip-placement="mouse" tooltip-animation="false" tooltip-append-to-body="true" ng-transclude></span>',
        link: function(scope, element, attrs) {
            element.css('color', '#6b6');
        }
    }
});

Application.Dna.directive('ligateBluntmatch', function() {
    return {
        restrict: 'EA',
        replace: false,
        transclude:true,
        template: '<span tooltip="Note! Blunt ends will only yield this product 50% of the time (fragment direction is random)" tooltip-placement="mouse" tooltip-animation="false" tooltip-append-to-body="true" ng-transclude></span>',
        link: function(scope, element, attrs) {
            element.css('color', '#6b6');
        }
    }
});
