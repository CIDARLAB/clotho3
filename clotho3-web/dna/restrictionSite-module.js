'use strict';

Application.Dna.service('Digest', ['Clotho', 'DNA', '$filter', function(Clotho, DNA, $filter) {

    var enzymes = {
        "BsaI" : {
            "name" : "BamHI",
            "match" : "ccannnnntgg",
            "cut" : "ccan|nnn^ntgg",
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
            "name" : "BamHI",
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
            "cut" : "c^tcga|g",
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
            "cut" : "g^gatc|c",
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
            "cut" : "g^aatt|c",
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
            "cut" : "t^ctag|a",
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
            "cut" : "a^ctag|t",
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
            "cut" : "c|tgca^g",
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
            "cut" : "a^agct|t",
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
            "cut" : "cag|nnn^ctg",
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
     * @returns {RegExp}
     */
    var createRegex = function (match) {
        return new RegExp(DNA.undegenerize(match), 'ig');
    };


    /**
     * @description Determine if an enzyme matches a given sequence
     * @param {String} sequence Sequence to check
     * @param {String} match String to use as Regexp, e.g. enzyme.match
     * @returns {boolean} Returns true if match, false otherwise
     */
    var sitePresent = function (sequence, match) {
        return (createRegex(match)).test(sequence);
    };

    var identifySite = function (match) {
        
    };


    /**
     * @description Find matches for a given oligo in a sequence
     * @param sequence
     * @param {String} match String to use as Regexp, e.g. enzyme.match
     * @returns {object} Matches in form { <match index> : <matched value> }
     */
    var findMatches = function (sequence, match) {
        var matches = {},
            cur,
            reg = createRegex(match);

        while ((cur = reg.exec(sequence)) != null) {
            matches[cur.index] = match[0];
        }

        return matches;
    };


    /**
     * @description Finds start indices of a oligo's match in a sequence
     * @param {String} sequence
     * @param {String} match String to use as Regexp, e.g. enzyme.match
     * @returns {Array} List of sites (start of match). Empty if no matches
     */
    var findIndices = function (sequence, match) {
        return Object.keys(findMatches(sequence, match));
    };


    /**
     * @description Adds a restriction site to the end of the sequence (3') with a given gap
     * @param {String} sequence
     * @param {Object} enzyme
     * @param {Number} gap number of nucleotides as spacer. Default is 3.
     * @param {boolean} fivePrime site should be added to 5' end (beginning) of sequence
     */
    var addRestrictionSite = function (sequence, enzyme, gap, fivePrime) {
        gap = angular.isDefined(gap) ? gap : 3;

        var adding = DNA.randomSequence(gap) + enzyme.match + DNA.randomSequence(gap);
        return (!!fivePrime) ? adding + sequence : sequence + adding;
    };

    var swapoutSites = function (sequence, enzyme) {
        //todo - need to find ORF

    };

    /**
     * @description Convert a sequence (fragment) cut by a restriction enzyme to have a sticky end.
     * @example given acgtacgACAGTG -> acgtacgA|CAGTG
     * @param {string} fragment
     * @param {object} enzyme
     * DEPRECATED
     */
    var addStickyEnd = function (fragment, enzyme) {
        fragment.replace(createRegex(enzyme.match), enzyme.cut);
    };



    var markSites = function (sequence, enzyme) {
        //add cut marks ^ and |
    };

    //purpose?
    var extractSites = function (fragment) {
        var findSiteReg = /[\^|\|].+?[\|\^]/gi,
            match,
            matches = {};

        while ((match = findSiteReg.exec(fragment)) != null) {
            matches[findSiteReg.lastIndex] = match;
        }
    };

    var makeCuts = function (sequence, enzyme) {
        //todo - handle enzymes like BsmBI with cut in form NNNNNNN (##/##)
        //todo - implement elsewhere

        //todo - handle multiple cuts (regExp.exec in loop)

        var indices = findIndices(sequence, enzyme.match),
            localReg = /(.+)^(.+)/ig,
            localCuts = localReg.exec(enzyme.cut),
            nonLocalReg = /\((\d+)\/(\d+)\)/ig,
            nonLocalCuts = nonLocalReg.exec(enzyme.cut),
            starts = [], ends = [];

        if (nonLocalCuts) {
            angular.forEach(indices, function (ind) {
                starts.push(ind + enzyme.match.length + nonLocalCuts[1]);
                ends.push(ind + enzyme.match.length + nonLocalCuts[2]);
            });

            //todo - finish, decide format for cuts

            return sequence;
        } else {
            return sequence.replace(createRegex(enzyme.match), enzyme.cut)
        }
    };





    /**
     * @description Determine fragments for a sequence cut by a single enzyme
     * @param {string} sequence
     * @param {object} enzyme
     * @returns {Array} Array of strings representing cut fragments
     */
    var simpleFragments = function (sequence, enzyme) {

        //todo - handle circular (beyond lookahead, should wrap around and join fragments)

        return sequence.split(createRegex(enzyme.match));
    };



    /* todo - rewrite so calculate:
     - if circular or linear
     - INFO NEED
        - if sites overlap - makes much harder to calculate
        - fragments
            - number
            - length of sequence
            - sticky ends (which enzyme)

     */
    var findCuts = function(sequence, enzymes) {

        //todo - finish rewriting
        var matches = {};

        angular.forEach(enzymes, function (enzyme) {

        });

        /**
         * @description Looks for sites only in one direction. Adds to locations array (extending object when revcomping may lead to overwrites)
         * @param {string} sequence
         * @param {Array} enzymes
         */
        function findSites (sequence, enzymes) {
            angular.forEach(enzymes, function (enz, index) {

                var seq = extendSequence(sequence, enz.match.length - 1);

                for (var index, offset = 0, search = angular.lowercase(seq);
                     (index = search.indexOf(angular.lowercase(enz.match), offset)) > -1;
                     offset = index + enz.match.length
                    ) {
                    //start
                    var start = index;
                    locations[start] ?
                        locations[start]['start'].push(enz) :
                        locations[start] = {"start" : [enz], "end" : []};
                    //end
                    var end = start + enz.match.length;
                    locations[end] ?
                        locations[end]['end'].push(enz) :
                        locations[end] = {"start" : [], "end" : [enz]};

                }
            });
        }


        findSites(sequence, enzymes);
        var revcomp = DNA.revcomp(sequence);
        findSites(revcomp, enzymes);
    
        console.log(locations);

        return locations

    };


    var getFragments = function (sequence, enzymes) {

        var locations = findCuts(sequence, enzymes),
            fragments = [],
            prevCut = undefined,
            revIndices = Object.keys(locations).reverse();

        //start from end, taking fragments
        angular.forEach(revIndices, function(index, ind) {
            angular.forEach(locations[index], function(array, key) {
                                        //this shouldn't happen
                if ((!array.length) || (sequence.length < index)) return;

                var seq = sequence;

                angular.forEach(array, function(enz) {
                    if (key == 'start') {

                        var frag = "";

                        //todo


                        //end sticky end
                        if (prevCut) {

                        }

                        frag = frag + sequence.substring(index);


                        //start sticky end
                        if (ind != revIndices.length - 1) {

                            //todo handle circular

                        }


                        fragments.unshift(frag);

                    }
                });
            });
        });

    };


    /**
     * @note Should only be called after fragments have been circularized
     * @note currently assumes field 'size' in fragment object, that fragments is an array of fragment objects
     */
    var sortFragments = function (fragments) {
        $filter('orderBy')(fragments, 'size');
    };


    return {
        enzymes : enzymes,

        extendSequence : extendSequence,
        circularize : circularize,
        createRegex : createRegex,
        sitePresent : sitePresent,
        findMatches : findMatches,
        findIndices : findIndices,

        addRestrictionSite : addRestrictionSite,
        addStickyEnd : addStickyEnd,
        makeCuts : makeCuts,
        extractSites : extractSites,
        simpleFragments : simpleFragments,

        findCuts: findCuts,             //in progress
        getFragments : getFragments,    //in progress

        sortFragments: sortFragments
    }

}]);

Application.Dna.directive('digestHighlight', ['Digest', '$filter', function(Digest, $filter) {
    return {
        restrict: 'A',
        require: 'ngModel',
        link : function(scope, element, attrs, ngModel) {

            var highlightSites = function (input) {
                //return $filter('highlight')(input, Digest.enzymes.BamHI.match, 'text-error');
                return Digest.makeCut(input, Digest.enzymes.BamHI)
            };

            var unhighlightSites = function (input) {
                return input.replace(/[^A-Z]/ig, '');
            };

            ngModel.$parsers.unshift(unhighlightSites);
            ngModel.$formatters.push(highlightSites);
        }
    }
}]);