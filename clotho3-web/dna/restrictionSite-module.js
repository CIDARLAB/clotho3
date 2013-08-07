'use strict';

Application.Dna.service('Digest', ['Clotho', 'DNA', '$filter', function(Clotho, DNA, $filter) {

    var enzymes = {
        "BsaI" : {
            "name" : "BasI",
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

    //todo - rewrite to handle new cut marks
    var regexps = {};
    regexps.localCuts = /\((.*?)([\^|\|])(.+?)([\|\^])(.*?)\)/ig;
    regexps.nonLocalCuts = /\((\d+)\/(\d+)\)/ig;
    regexps.findOverhang = /([\^|\|])(.+?)([\|\^])/ig;

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
        return (createRegex(match, true)).test(sequence);
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


    var swapoutSites = function (sequence, enzyme, orfOffset) {
        var offset = orfOffset ? orfOffset : DNA.probableORF(sequence);

        // todo
    };


    /**
     * @description Find matches for a given oligo in a sequence
     * @param sequence
     * @param {String} match String to use as Regexp, e.g. enzyme.match
     * @param {boolean} bothStrands Whether revcomp of match should be searched as well
     * @returns {object} Matches in form { <match index> : <matched value> }
     */
    var findMatches = function (sequence, match, bothStrands) {
        var matches = {},
            cur,
            reg = createRegex(match, bothStrands !== false);

        while ((cur = reg.exec(sequence)) != null) {
            matches[cur.index] = match[0];
        }

        return matches;
    };


    /**
     * @description Finds start indices of a oligo's match in a sequence
     * @param {String} sequence
     * @param {String} match String to use as Regexp, e.g. enzyme.match
     * @param {boolean} bothStrands Whether revcomp of match should be searched as well
     * @returns {Array} List of sites (start of match). Empty if no matches
     */
    var findIndices = function (sequence, match, bothStrands) {
        return Object.keys(findMatches(sequence, match, bothStrands));
    };


    /**
     * @description Marks enzyme recognitions sites and cut marks on a sequence
     * @param {string} sequence
     * @param {object} enzyme
     * @returns {string} Returns sequence with cuts and recognition sites marked
     *
     * @example BamHI { "cut" : "g^gatc|c" } NNNNNggatccNNNN -> NNNNN(g^gatc|c)NNNNN
     * @example BsmbI { "cut" : "cgtctc (1/5)" } NNNcgtctcNNNNNNNNN -> NNN(cgtctc)N^NNNN|NNNN
     * @example BsmbI { "cut" : "cgtctc (1/5)" } NNNNNNNNNgagacgNNN -> NNNN|NNNN^N(gagacg)NNN
     *
     * note - won't handle case where enzyme cuts on either side of reg site (i.e. one negative, one positive number
     */
    var markSites = function (sequence, enzyme) {

        //todo - match cut (i.e. ^ vs. | being first)

        //todo - handle 3' overhangs (e.g. ATCGTG (20/18))
        //todo - handle one positive one negative number

        sequence = DNA.dnaOnly(sequence);

        var nonLocalCuts = (/\((\d+)\/(\d+)\)/ig).exec(enzyme.cut);

        if (nonLocalCuts) {

            //matches to extract (for constructing regexp)
            var brs = {};
            brs.enz = '(' + enzyme.match + ')';
            brs.rev = '(' + DNA.revcomp(enzyme.match) + ')';


            //cuts before recognition sequence
            if (nonLocalCuts[1] < 0) {
                brs.gap1 = '(.{'+Math.abs(nonLocalCuts[1]-nonLocalCuts[2])+'})';
                brs.gap2 = '(.{'+Math.abs(nonLocalCuts[2])+'})';

                //forward
                var cutFor = new RegExp(brs.gap2 + brs.gap1 + brs.enz, 'ig');
                sequence = sequence.replace(cutFor, function(match, $1, $2, $3, off, orig) {
                    return [ '|' + $1 + '^' + $2 + '(' + $3 + ')']
                });
                //reverse
                var cutRev = new RegExp(brs.enz + brs.gap1 + brs.gap2, 'ig');
                sequence = sequence.replace(cutRev, function(match, $1, $2, $3, off, orig) {
                    return [ '(' + $1 + ')' + $2 + '^' + $3 + '|']
                });

            } else {
                brs.gap1 = '(.{'+nonLocalCuts[1]+'})';
                brs.gap2 = '(.{'+(nonLocalCuts[2]-nonLocalCuts[1])+'})';

                //forward
                var cutFor = new RegExp(brs.enz + brs.gap1 + brs.gap2, 'ig');
                sequence = sequence.replace(cutFor, function(match, $1, $2, $3, off, orig) {
                    return [ '(' + $1 + ')' + $2 + '^' + $3 + '|']
                });
                //reverse
                var cutRev = new RegExp(brs.gap2 + brs.gap1 + brs.rev, 'ig');
                sequence = sequence.replace(cutRev, function(match, $1, $2, $3, off, orig) {
                    return [ '|' + $1 + '^' + $2 + '(' + $3 + ')']
                });
                
                console.log(sequence);
            }

        } else {
            //need to check for forward and reverse matches, update replacement appropriately
            sequence = sequence.replace(createRegex(enzyme.match), '(' + enzyme.cut + ')');
            sequence = sequence.replace(createRegex(DNA.revcomp(enzyme.match)), '(' + DNA.revcomp(enzyme.cut) + ')');
        }

        return sequence;
    };

    /**
     * @description Marks recognition sites by surrounding with parentheses
     * @param {string} sequence
     * @param {string} enzyme
     * @returns {string} Sequence with marked restriction sites
     */
    var markMatches = function (sequence, enzyme) {
        return sequence.replace(createRegex(enzyme.match, true), function(match) {
            return '(' + match + ')'
        });
    };


    var markCuts = function (sequence, enzyme) {
        sequence = markSites(sequence, enzyme);
        return sequence.replace(/\(|\)/gi, '');
    };


    var extractSites = function (sequence) {
        var reg = regexps.findOverhang,
            match,
            matches = {};

        while ((match = (reg).exec(sequence)) != null) {
            matches[reg.lastIndex] = match;
        }

        return matches;
    };

    /**
     * @description Determine fragments for a sequence cut by a single enzyme
     * @param {string} sequence
     * @param {object} enzyme
     * @param {boolean=} Whether fragments should be circularized. Defaults to true
     * @returns {Array} Array of strings representing cut fragments
     */
    var makeCuts = function (sequence, enzyme, circular) {
        //todo - more complex logic?
        var fragments = (markSites(sequence, enzyme)).split('^');
        console.log(fragments);

        return (circular === false) ? fragments : circularize(fragments);
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
        identifySite : identifySite,
        addRestrictionSite : addRestrictionSite,
        swapoutSites : swapoutSites,

        findMatches : findMatches,
        findIndices : findIndices,
        markSites : markSites,
        markMatches : markMatches,
        markCuts : markCuts,
        extractSites : extractSites,
        makeCuts : makeCuts,

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
                return Digest.markSites(input, Digest.enzymes.BsmbI)
            };

            var unhighlightSites = function (input) {
                return input.replace(/[^A-Z]/ig, '');
            };

            ngModel.$parsers.unshift(unhighlightSites);
            ngModel.$formatters.push(highlightSites);
        }
    }
}]);