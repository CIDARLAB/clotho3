'use strict';

Application.Dna.service('Digest', ['Clotho', 'DNA', '$filter', function(Clotho, DNA, $filter) {

    var enzymes = {
        "BsaI" : {
            "name" : "BasI",
            "match" : "ccannnnntgg",
            "cut" : "ccan_nnn^ntgg",
            "strand" : "",
            "methylation" : false,
            "overhang" : 3,
            "type" : "II",
            "subtype" : "P",
            "notes" : {},
            "buffer" : "",
            "star activity" : false,
            "comment" : "",
            "rebase" : "http://rebase.neb.com/rebase/enz/BsaI.html",
            "personal" : {},
            "citations" : {},
            "ordering" : {}
        },
        "BsmbI" : {
            "name" : "BsmbI",
            "match" : "cgtctc",
            "cut" : "cgtctc (1/5)",
            "strand" : "",
            "methylation" : true,
            "overhang" : 4,
            "type" : "II",
            "subtype" : "S",
            "notes" : {},
            "buffer" : "",
            "star activity" : false,
            "comment" : "",
            "rebase" : "http://rebase.neb.com/rebase/enz/BsmbI.html",
            "personal" : {},
            "citations" : {},
            "ordering" : {}
        },
        "XhoI" : {
            "name" : "XhoI",
            "match" : "ctcgag",
            "cut" : "c^tcga_g",
            "strand" : "",
            "methylation" : true,
            "overhang" : 4,
            "type" : "II",
            "subtype" : "P",
            "notes" : {},
            "buffer" : "",
            "star activity" : false,
            "comment" : "flanking C5 methylation can slow cleavage, ATCTCTCGAGTCTA is cut v. slowly",
            "rebase" : "http://rebase.neb.com/rebase/enz/XhoI.html",
            "personal" : {},
            "citations" : {},
            "ordering" : {}
        },
        "BamHI" : {
            "name" : "BamHI",
            "match" : "ggatcc",
            "cut" : "g^gatc_c",
            "strand" : "",
            "methylation" : {},
            "overhang" : 4,
            "type" : "II",
            "subtype" : "P",
            "notes" : {},
            "buffer" : "",
            "star activity" : false,
            "comment" : "",
            "rebase" : "http://rebase.neb.com/rebase/enz/BamHI.html",
            "personal" : {},
            "citations" : {},
            "ordering" : {}
        },
        "EcoRI" : {
            "name" : "EcoRI",
            "match" : "gaattc",
            "cut" : "g^aatt_c",
            "strand" : "",
            "methylation" : {},
            "overhang" : 4,
            "type" : "II",
            "subtype" : "P",
            "notes" : {},
            "buffer" : "",
            "star activity" : true,
            "comment" : "",
            "rebase" : "http://rebase.neb.com/rebase/enz/EcoRI.html",
            "personal" : {},
            "citations" : {},
            "ordering" : {}
        },
        "XbaI" : {
            "name" : "XbaI",
            "match" : "tctaga",
            "cut" : "t^ctag_a",
            "strand" : "",
            "methylation" : {},
            "overhang" : 4,
            "type" : "II",
            "subtype" : "P",
            "notes" : {},
            "buffer" : "",
            "star activity" : true,
            "comment" : "",
            "rebase" : "http://rebase.neb.com/rebase/enz/XbaI.html",
            "personal" : {},
            "citations" : {},
            "ordering" : {}
        },
        "SpeI" : {
            "name" : "XbaI",
            "match" : "actagt",
            "cut" : "a^ctag_t",
            "strand" : "",
            "methylation" : true,
            "overhang" : 4,
            "type" : "II",
            "subtype" : "P",
            "notes" : {},
            "buffer" : "",
            "star activity" : true,
            "comment" : "",
            "rebase" : "http://rebase.neb.com/rebase/enz/SpeI.html",
            "personal" : {},
            "citations" : {},
            "ordering" : {}
        },
        "PstI" : {
            "name" : "PstI",
            "match" : "ctgcag",
            "cut" : "c_tgca^g",
            "strand" : "",
            "methylation" : {},
            "overhang" : 4,
            "type" : "II",
            "subtype" : "P",
            "notes" : {},
            "buffer" : "",
            "star activity" : true,
            "comment" : "",
            "rebase" : "http://rebase.neb.com/rebase/enz/PstI.html",
            "personal" : {},
            "citations" : {},
            "ordering" : {}
        },
        "HindIII" : {
            "name" : "HindIII",
            "match" : "aagctt",
            "cut" : "a^agct_t",
            "strand" : "",
            "methylation" : {},
            "overhang" : 4,
            "type" : "II",
            "subtype" : "P",
            "notes" : {},
            "buffer" : "",
            "star activity" : true,
            "comment" : "",
            "rebase" : "http://rebase.neb.com/rebase/enz/HindIII.html",
            "personal" : {},
            "citations" : {},
            "ordering" : {}
        },
        "AlwnI" : {
            "name" : "AlwnI",
            "match" : "cagnnnctg",
            "cut" : "cag_nnn^ctg",
            "strand" : "",
            "methylation" : {},
            "overhang" : 3,
            "type" : "II",
            "subtype" : "P",
            "notes" : {},
            "buffer" : "",
            "star activity" : false,
            "comment" : "",
            "rebase" : "http://rebase.neb.com/rebase/enz/AlwnI.html",
            "personal" : {},
            "citations" : {},
            "ordering" : {}
        },
        "AcuI" : {
            "name" : "AcuI",
            "match" : "ctgaag",
            "cut" : "ctgaag (16/14)",
            "strand" : "",
            "methylation" : {},
            "overhang" : 2,
            "type" : "II",
            "subtype" : "G, S, alpha",
            "notes" : {},
            "buffer" : "",
            "star activity" : false,
            "comment" : "Methylated product: 6-methyladenosine (base undetermined)",
            "rebase" : "http://rebase.neb.com/rebase/enz/XbaI.html",
            "personal" : {},
            "citations" : {},
            "ordering" : {}
        },
        "BseRI" : {
            "name" : "BseRI",
            "match" : "gaggag",
            "cut" : "gaggag (10/8)",
            "strand" : "",
            "methylation" : {},
            "overhang" : 2,
            "type" : "II",
            "subtype" : "G, S",
            "notes" : {},
            "buffer" : "",
            "star activity" : false,
            "comment" : "Methylated Product: 6-methyladenosine (base undetermined)",
            "rebase" : "http://rebase.neb.com/rebase/enz/XbaI.html",
            "personal" : {},
            "citations" : {},
            "ordering" : {}
        }
    };

    //calculated on the fly, reverse lookup by match or cut sequence
    var enzymesReverse = {};
    angular.forEach(enzymes, function (enz, name) {
        enzymesReverse[enz.match] = name;
        enzymesReverse[enz.cut] = name;
    });

    var cuts = {
        'blunt' : '|',
        'main'  : '^',
        'comp'  : '_'
    };

    var regexps = {};
    regexps.localCuts = /\(((.*?)([\^|_])(.+?)([\^|_])(.*?))|((.*?)(\|)(.*?))\)/ig;
    regexps.nonLocalCuts = /\((\d+)\/(\d+)\)/ig;
    regexps.findCut = /(\|)|([\^_])(.+?)([_\^])/ig;
    regexps.findBlunt = /(\|)/ig;
    regexps.findOverhang = /([\^_])(.+?)([_\^])/ig;

    /**
     * @description Extend a sequence by a defined length, repeating monomers from the beginning.
     * @param {String} sequence
     * @param {Number} lookahead Length to extend. Example lookahead length is (enzyme.match.length - 1)
     * @returns {string} Extended sequence
     */
    var extendSequence = function (sequence, lookahead) {
        return (sequence + sequence.substr(lookahead));
    };

    /**
     * @description Given an UNSORTED (i.e. order in which linear sequence was cut) array of fragments, will join the first fragment to the last fragment, approximating the result of cutting a circular instead of linear sequence.
     * @param {Array} fragments Array of sequence fragments, MUST BE ORDERED SUCH THAT fragments[0] is the first one, and fragments[length-1] is the last one, if the sequence were linear.
     * @returns {Array}
     */
    var circularize = function (fragments) {
        if (fragments.length > 1) {
            var first = fragments.shift();
            fragments[fragments.length-1] = fragments[fragments.length-1] + first;
        }
        return fragments;
    };


    /**
     * @description Creates a regular expression for a given enzyme, converting degenerate match sequence appropriately
     * @example ACGTBDHVN -> ACGT[CGTU][AGTU][ACTU][ACG][ACGTU]
     * @param {String} match String to use as Regexp, e.g. enzyme.match
     * @param {boolean=} doubleStrand  Whether match should account for both strands (revcomp). Default is false.
     * @returns {RegExp}
     */
    var createRegex = function (match, doubleStrand) {
        var degen = doubleStrand ? (match + '|' + DNA.revcomp(match)) : match;
        return new RegExp(DNA.undegenerize(degen), 'ig');
    };


    /**
     * @description Determine if an enzyme matches a given sequence
     * @param {String} sequence Sequence to check
     * @param {String} match String to use as Regexp, e.g. enzyme.match
     * @returns {boolean} Returns true if match, false otherwise
     */
    var sitePresent = function (sequence, match) {
        return (createRegex(match, true)).test(extendSequence(sequence, match.length-1));
    };


    /**
     * @description Identify the enzyme whose recognition sequence matches 'match'
     * @param match Should match either enzyme.match or enzyme.cut
     * @returns {string} Name of enzyme, same as key in enzymes object
     */
    var identifySite = function (match) {
        return enzymesReverse[match];
    };

    /**
     * @description Adds a restriction site to the end of the sequence (3') with a given gap
     * @param {String} sequence
     * @param {Object} enzyme
     * @param {Number} gap Number of nucleotides as spacer. Default is 3.
     * @param {boolean} fivePrime Whether site should be added to 5' (beginning) or 3' (end)
     * @param {boolean} oppStrand Whether site is on opposite strand (should be reverse complemented)
     */
    var addRestrictionSite = function (sequence, enzyme, gap, fivePrime, oppStrand) {
        gap = angular.isDefined(gap) ? gap : 3;

        var adding = DNA.randomSequence(gap) +
            (oppStrand ? DNA.revcomp(enzyme.match) : enzyme.match) +
            DNA.randomSequence(gap);
        return (!!fivePrime) ? adding + sequence : sequence + adding;
    };

    /**
     * @description attempts to alter sequence preserving protein sequence, but removing sequence specified (e.g. enzyme.match)
     * @param sequence
     * @param {RegExp|String} match Either a Regular Expression, or string to use as Regexp (e.g. enzyme.match)
     * @param forceOffset
     * @return {boolean|string} Altered sequence, false if could not swap. Changes will be lowercase
     */
    var swapoutSites = function (sequence, match, forceOffset) {
        var offset = forceOffset ? forceOffset : DNA.probableORF(sequence),
            matches = findMatches(sequence, match, true),
            indices = Object.keys(matches);

        //we're ok
        if (matches.length == 0) return sequence;

        var testReg = (match instanceof RegExp) ? match : createRegex(match, true);

        //loop each index with a match
        loop_eachIndex: for (var i = 0; i < indices.length; i++) {
            var index = indices[i],
                curMatch = matches[index],
                start = index - (index % 3) + offset,
                numCodons = Math.ceil((curMatch.length + (index - start)) % 3);

            //loop through codons containing match
            loop_eachCodon : for (var j = 0; j <= numCodons; j++) {
                var curIndex = start + j*3,
                    codon = sequence.substr(curIndex, 3).toLowerCase(),
                    amino = DNA.complements.translate[codon],
                    options = DNA.complements.reverse_translate[amino],
                    replacement;

                if (options.length <= 1)
                    continue;

                // loop through options for codon.
                // If remove match, break loop_eachCodon. Otherwise, return false;
                loop_eachOption : for (var k = 0; k < options.length; k++) {
                    if (options[k] == codon) continue;



                    //todo - don't test whole sequence, just this selection
                    var testSeq = sequence.substr(0, curIndex) + options[k] + sequence.substr(curIndex+3);
                    if (!testReg.test(testSeq)) {
                        sequence = testSeq;
                        break loop_eachCodon;
                    }

                }



                //todo - finish -- if not good after each codon, return false
                if (testReg.test()) {}






                //future - check codon frequencies


            }

            // todo
            //handle degenerate matches


            /*
            offset 1
            acgtacgtATGCATacgtacgt

            index 8

            -> codons @ 7-8-9 (tAT), 10-11-12 (GCA), 13-14-15 (Tac)
                start = index - (index % 3) + offset = 7
                numCodons = Math.ceil((match.length + (index - start)) % 3)

            ->



             */


            if (impossible) return false
        }

        return sequence
    };


    /**
     * @description Find matches for a given oligo in a sequence
     * @param sequence
     * @param {RegExp|String} match Either a Regular Expression, or string to use as Regexp (e.g. enzyme.match)
     * @param {boolean} bothStrands Whether revcomp of match should be searched as well. Default is true
     * @returns {object} Matches in form { <match index> : <matched value> }
     */
    var findMatches = function (sequence, match, bothStrands) {
        var matches = {},
            cur,
            reg = (match instanceof RegExp) ? match : createRegex(match, bothStrands !== false);

        //future - move to fuzzy search

        while ((cur = reg.exec(extendSequence(sequence, match.length-1))) != null) {
            console.log(match);
            matches[cur.index] = match;
        }

        return matches;
    };


    /**
     * @description Finds start indices of a oligo's match in a sequence
     * @param {String} sequence
     * @param {RegExp|String} match Either a Regular Expression, or string to use as Regexp (e.g. enzyme.match)
     * @param {boolean} bothStrands Whether revcomp of match should be searched as well. Default is true
     * @returns {Array} List of sites (start of match). Empty if no matches
     */
    var findIndices = function (sequence, match, bothStrands) {
        return Object.keys(findMatches(sequence, match, bothStrands));
    };



    var removeCuts = function(sequence) {
        return sequence.replace(/[\|\^_]/gi, '')
    };

    var removeMatches = function(sequence) {
        return sequence.replace(/[\(\)]/gi, '')
    };

    var removeMarks = function (sequence) {
        return removeMatches(removeCuts(sequence));
    };

    /**
     * @description Marks enzyme recognitions sites and cut marks on a sequence
     * @param {string} sequence
     * @param {object} enzyme
     * @returns {string} Returns sequence with cuts and recognition sites marked
     *
     * @example [Normal]
     * BamHI { "cut" : "g^gatc_c" }
     * NNNNNggatccNNNN -> NNNNN(g^gatc_c)NNNNN
     *
     * @example [Nonlocal]
     * BsmbI { "cut" : "cgtctc (1/5)" }
     * NNNcgtctcNNNNNNNNN -> NNN(cgtctc)N^NNNN_NNNN
     *
     * @example [Nonlocal Revcomp]
     * BsmbI { "cut" : "cgtctc (1/5)" }
     * NNNNNNNNNgagacgNNN -> NNNN_NNNN^N(gagacg)NNN
     *
     * @example [3' overhang]
     * BseRI { "cut" : "gaggag (10/8)" }
     * NNNgaggagNNNNNNNNNNNNNNNNNNNN -> NNN(gaggag)NNNNNNNN_NN^NNNNNNNNNN
     *
     * @example [3' side]
     * {"cut" : "atgcat (-5/-1)" }
     * NNNNNNNNNNatgcatNNN -> NNNNN^NNNN_N(atgcat)NNN
     *
     * @example [3' side, revcomp]
     * {"cut" : "atgcat (-5/-1)" }
     * NNNatgcatNNNNNNNNNN -> NNN(atgcat)N^NNNN_NNNNN
     */

        //future - handle cuts on both sides (e.g. Bsp24I (8/13)GACNNNNNNTGG(12/7) and one on either side

        //todo - multiple enzymes
    var markSites = function (sequence, enzyme) {

        if (!enzyme) return sequence;

        sequence = DNA.dnaOnly(sequence);

        var nonLocalCuts = (/\((\d+)\/(\d+)\)/ig).exec(enzyme.cut);

        if (nonLocalCuts) {

            //matches to extract (for constructing regexp) - backreferences
            var brs = {};
            brs.enz = '(' + enzyme.match + ')';
            brs.rev = '(' + DNA.revcomp(enzyme.match) + ')';

            var cut53 = ['^', '_'],
                cut35 = ['_', '^'];

            // whether first cut is parent or complimentary strand (^..._ vs _...^)
            // i.e. 3' or 5' overhang
            var cut = (nonLocalCuts[1] < nonLocalCuts[2]) ? cut53 : cut35;

            //cuts before recognition sequence (e.g. AGTACT (-5/-1)
            if (nonLocalCuts[1] < 0) {

                brs.gap1 = '(.{'+Math.abs(nonLocalCuts[1]-nonLocalCuts[2])+'})';
                brs.gap2 = '(.{'+Math.abs(nonLocalCuts[2])+'})';

                //forward
                var cutFor = new RegExp(brs.gap2 + brs.gap1 + brs.enz, 'ig');
                sequence = sequence.replace(cutFor, function(match, $1, $2, $3, off, orig) {
                    return [ cut[0] + $1 + cut[1] + $2 + '(' + $3 + ')']
                });
                //reverse
                var cutRev = new RegExp(brs.enz + brs.gap1 + brs.gap2, 'ig');
                sequence = sequence.replace(cutRev, function(match, $1, $2, $3, off, orig) {
                    return [ '(' + $1 + ')' + $2 + cut[1] + $3 + cut[0]]
                });

            }
            else {
                brs.gap1 = '(.{'+nonLocalCuts[1]+'})';
                brs.gap2 = '(.{'+(nonLocalCuts[2]-nonLocalCuts[1])+'})';

                //forward
                var cutFor = new RegExp(brs.enz + brs.gap1 + brs.gap2, 'ig');
                sequence = sequence.replace(cutFor, function(match, $1, $2, $3, off, orig) {
                    return [ '(' + $1 + ')' + $2 + cut[0] + $3 + cut[1]]
                });
                //reverse
                var cutRev = new RegExp(brs.gap2 + brs.gap1 + brs.rev, 'ig');
                sequence = sequence.replace(cutRev, function(match, $1, $2, $3, off, orig) {
                    return [ cut[1] + $1 + cut[0] + $2 + '(' + $3 + ')']
                });
                
                console.log(sequence);
            }

        } else {
            //note - already account for cut marks being reverse complimented
            sequence = sequence.replace(createRegex(enzyme.match), '(' + enzyme.cut + ')');
            sequence = sequence.replace(createRegex(DNA.revcomp(enzyme.match)), '(' + DNA.revcomp(enzyme.cut) + ')');
        }


        //handle lookahead
        var wrap = sequence.substring(sequence.length - enzyme.match.length) + sequence.substr(enzyme.match.length-1);
        if (createRegex(enzyme.match).exec(wrap)) {

            //todo - handle, including nonLocal

            if (nonLocalCuts) {

            } else {

            }
        }

        return sequence;
    };

    /**
     * @description Marks recognition sites by surrounding with parentheses
     * @param {string} sequence
     * @param {object} enzyme
     * @returns {string} Sequence with marked restriction sites
     */
    var markMatches = function (sequence, enzyme) {
        return removeCuts(markSites(sequence, enzyme))
    };


    /**
     * @description Marks only cut sites without marking matched sequences
     * @param {string} sequence
     * @param {object} enzyme
     * @returns {string}
     */
    var markCuts = function (sequence, enzyme) {
        return removeMatches(markSites(sequence, enzyme))
    };

    /**
     * @description
     * @param sequence
     * @returns {object} key is index of match, value is array...
     * Some keys:
     *      0-4 : <matches from regex>
     *      index: index of match
     *      input: passed in sequence
     *      isBlunt: if blunt cut (i.e. "|")
     *      length: length of match
     *      match: matched overhang sequence, including marks
     *      rev : true for 3' overhang (e.g. nnn_nnnn^nnn)
     *
     * @example "acgt^ct_acagctagcta|gctagctagct_cgta^agagctacga"
     4: Array[5]
         0: "^ct_"
         1: undefined
         2: "^"
         3: "ct"
         4: "_"
         index: 4
         input: "acgt^ct_acagctagcta|gctagctagct_cgta^agagctacga"
         isBlunt: false
         length: 4
         match: "^ct_"
         rev: false
     19: Array[5]
         0: "|"
         1: "|"
         2: undefined
         3: undefined
         4: undefined
         index: 19
         input: "acgt^ct_acagctagcta|gctagctagct_cgta^agagctacga"
         isBlunt: true
         length: 1
         match: "|"
         rev: null
     31: Array[5]
         0: "_cgta^"
         1: undefined
         2: "_"
         3: "cgta"
         4: "^"
         index: 31
         input: "acgt^ct_acagctagcta|gctagctagct_cgta^agagctacga"
         isBlunt: false
         length: 6
         match: "_cgta^"
         rev: true
     */
    var findOverhangs = function (sequence) {
        var regCut = regexps.findCut,
            match,
            matches = {};

        //testing
        //sequence = "acgt^ct_acagctagcta|gctagctagct_cgta^agagctacga";

        while ((match = (regCut).exec(extendSequence(sequence, match.length-1))) != null) {
            match.match = match[0];
            match.isBlunt = (match.match == '|');
            match.rev = match.isBlunt ? null : (match[2] == '^');
            match.length = match.match.length;
            matches[match.index] = match;
        }

        return matches;
    };

    /**
     * @description Determine fragments for a sequence cut by a single enzyme
     * @param {string} sequence
     * @param {object=} enzymes If none passed, assume sequence is marked
     * @param {boolean=} circular Whether fragments should be circularized. Default: true
     * @returns {Array} Array of strings representing cut fragments
     */

        /*
        nnnn^nn_nnnn -> [nnnn^nn_, ^nn_nnnn]
        nnn|nnn -> [nnn, nnn]
         */
    var makeCuts = function (sequence, enzymes, circular) {

        if (!!enzymes)
            sequence = markCuts(sequence, enzymes);

        var lastIndex = 0,
            lastMark = "",
            keys = [],
            fragments = [],
            overhangs = findOverhangs(sequence);

        //ensure sorted
        for (var ind in overhangs) {
            if (overhangs.hasOwnProperty(ind)) {
                keys.push(ind)
            }
        }
        keys.sort();

        //split into fragments
        for (var i = 0; i < keys.length; i++) {
            var mark = overhangs[keys[i]],
                newMark = (mark.isBlunt) ? '' : mark.match,
                frag = lastMark + sequence.substring(lastIndex, mark.index) + newMark;

            fragments.push(frag);

            lastIndex = mark.index;
            lastMark = newMark;
        }

        return (circular === false) ? fragments : circularize(fragments);
    };


    return {
        //config
        enzymes : enzymes,

        //utility
        extendSequence : extendSequence,
        createRegex : createRegex,
        circularize : circularize,

        //site basics
        sitePresent : sitePresent,
        identifySite : identifySite,
        findMatches : findMatches,
        findIndices : findIndices,

        //find and mark
        markSites : markSites,
        markMatches : markMatches,
        markCuts : markCuts,
        removeCuts : removeCuts,
        removeMatches : removeMatches,
        removeMarks : removeMarks,
        findOverhangs : findOverhangs,

        //manipulation
        makeCuts : makeCuts,
        addRestrictionSite : addRestrictionSite,
        swapoutSites : swapoutSites

    }
}]);

Application.Dna.directive('digestMark', ['Digest', '$filter', '$parse', function(Digest, $filter, $parse) {
    return {
        restrict: 'A',
        require: 'ngModel',
        link : function(scope, element, attrs, ngModel) {

            var enz;
            scope.$watch(attrs['digestMark'], function (val) {
                enz = val;
                ngModel.$render(); //fixme doens't work
            });

            var highlightSites = function (input) {
                //return $filter('highlight')(input, Digest.enzymes.BamHI.match, 'text-error');
                return Digest.markSites(input, enz)
            };

            var unhighlightSites = function (input) {
                return Digest.removeMarks(input);
            };

            ngModel.$parsers.unshift(unhighlightSites);
            ngModel.$formatters.push(highlightSites);
        }
    }
}]);

//todo - migrate Plasmid directive into this module
Application.Dna.directive('digestHighlight', ['Digest', '$parse', '$compile', function(Digest, $parse, $compile) {
    return {
        restrict: 'A',
        require: 'ngModel',

        link : function(scope, element, attrs, ngModel) {

            scope.$watch(attrs['digestHighlight'], function (val) {
                console.log(val);
                scope.highlightEnz = val;
                reprocess();
            });

            scope.$watch(function() {
                return ngModel.$modelValue
            }, function(newval, oldval) {
                reprocess()
            });
            
            function reprocess() {
                var seqSites = Digest.markSites(ngModel.$modelValue, scope.highlightEnz);
                var findMatch = /\((.+?)\)/gi;
                var addedAnnotations = seqSites.replace(findMatch, '<digest-annotation>$1</digest-annotation>');
                element.html($compile('<div>' + addedAnnotations + '</div>')(scope))
            }
            

            //todo - show cut marks nicer

        }
    }

}]);

Application.Dna.directive('digestAnnotation', ['$tooltip', function($tooltip) {

    return {
        restrict : 'EA',
        replace: false,
        transclude:true,
        template: '<span tooltip="{{ highlightEnz.name }}" tooltip-placement="mouse" tooltip-animation="false" tooltip-append-to-body="true" ng-transclude></span>',
        compile: function compile(tElement, tAttrs, transclude) {
        return {
            pre: function preLink(scope, element, attrs, controller) {

            },
            post: function(scope, element, attrs, ctrl) {
                console.log(scope.highlightEnz);

                element.css({backgroundColor: '#fcc'});

            }
        }
    }
}


}])
