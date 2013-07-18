'use strict';

//todo - functionality
// directives that transform don't affect ngModel on change
// codon frequency e.g. e-coli --->> smart backTranslate
// gibbs calculation
// restriction site service
// FASTA reader, genbank
// alignment
// mutation
// pcr products
// fuzzy search

Application.Dna.service('DNA', ['$filter', function($filter) {

    /*
     iupac_nucleotides =
     [ 'A', 'T', 'C', 'G', 'U',    Canonical bases
     'B', 'V', 'D', 'H',           B = ^A, V= ^T, D= ^C, H= ^G
     'S', 'W',                     Strong and Weak: GC vs. AT
     'K', 'M',                     "Keto" and "aMino": GT vs. AC
     'R', 'Y',                     "puRine" and "pYrimidine": AG vs. CT
     'N' ]                         Wildcards: any N
     */

    var regexps = {},
        monomers = {},
        maps = {},
        complements = {},
        frequencies = {},
        weights = {};

    regexps.letters = /[A-Z]/ig;
    regexps.nucleotide = /[gautcGAUTC]/g;
    regexps.dna = /[gatcGATC]/g;
    regexps.rna = /[gaucGAUC]/g;
    regexps.dna_degnerate = /[gatucryswkmbdhvnxGATUCRYSWKMBDHVNX]/g;
    regexps.protein = /[ACDEFGHIKLMNPQRSTVWYZacdefghiklmnpqrstvwyz\*]/g;
    regexps.protein_degenerate = /[ABCDEFGHIKLMNPQRSTVWYXZabcdefghiklmnpqrstvwyxz\*]/g;

    regexps.strip = {};
    regexps.strip.dna = /[^gatcGATC]/g;
    regexps.strip.dna_degenerate = /[^gatucryswkmbdhvnxGATUCRYSWKMBDHVNX]/g;
    regexps.strip.rna = /[^gaucGAUC]/g;
    regexps.strip.protein = /[^ACDEFGHIKLMNPQRSTVWYZacdefghiklmnpqrstvwyz\*]/g;
    regexps.strip.protein_degenerate = /[^ABCDEFGHIKLMNPQRSTVWYXZabcdefghiklmnpqrstvwyxz\*]/g;
    regexps.strip.alignment = /[^\.\-]/g;
    regexps.strip.non_letters = /[^A-Z]/ig;
    regexps.strip.whitespace = /\s/g;


    monomers.dna = ['a', 'c', 't', 'g'];
    monomers.rna = ['a', 'c', 'u', 'g'];
    monomers.nucleotide = ['a', 'c', 't', 'u', 'g'];
    monomers.dna_degenerate = ['g','a','t','u','c','r','y','s','w','k','m','b','d','h','v','n','x'];
    monomers.protein = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','Z','*'];
    monomers.protein_degnerate = ['A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z','*'];
    //dna, rna, protein
    monomers.all = ['A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','U','V','W','X','Y','Z','*'];


    maps.nucleotide_degenerate = {
        'A': 'A', 'B': '[CGTU]', 'C': 'C', 'D': '[AGTU]', 'G': 'G', 'H': '[ACTU]',
        'K': '[GTU]', 'M': '[AC]', 'N': '[ACGTU]', 'R': '[AG]', 'S': '[CG]',
        'T': 'T', 'U': 'U', 'V': '[ACG]', 'W': '[ATU]', 'Y': '[CTU]',
        '.': '[ACGTU]', '-': '[ACGTU]'};
    maps.dna_degenerate = {'A':'A','B':'[CGT]','C':'C','D':'[AGT]','G':'G',
        'H':'[ACT]','K':'[GT]','M':'[AC]','N':'[ACGT]','R':'[AG]','S':'[CG]',
        'T':'T','V':'[ACG]','W':'[AT]','Y':'[CT]','.':'[ACGT]','-':'[ACGT]'};
    maps.rna_degenerate = {'A':'A','B':'[CGU]','C':'C','D':'[AGU]','G':'G',
        'H':'[ACU]','K':'[GU]','M':'[AC]','N':'[ACGU]','R':'[AG]','S':'[CG]',
        'U':'U','V':'[ACG]','W':'[AU]','Y':'[CU]','.':'[ACGU]','-':'[ACGU]'};

    maps.amino_one_to_three = {
        'A': 'Ala',
        'B': 'Asx',
        'C': 'Cys',
        'D': 'Asp',
        'E': 'Glu',
        'F': 'Phe',
        'G': 'Gly',
        'H': 'His',
        'I': 'Ile',
        'K': 'Lys',
        'L': 'Leu',
        'M': 'Met',
        'N': 'Asp',
        'P': 'Pro',
        'Q': 'Gln',
        'R': 'Arg',
        'S': 'Ser',
        'T': 'Thr',
        'V': 'Val',
        'W': 'Trp',
        'X': 'X',
        'Y': 'Tyr',
        'Z': 'Glx',
        '*': '*'
    };
    maps.amino_one_to_full = {
        'A': 'Alanine',
        'B': '[Aspartic Acid or Asparagine]',
        'C': 'Cysteine',
        'D': 'Aspartic Acid',
        'E': 'Glutamic Acid',
        'F': 'Phenylalanine',
        'G': 'Glycine',
        'H': 'Histidine',
        'I': 'Isoleucine',
        'K': 'Lysine',
        'L': 'Leucine',
        'M': 'Methionine',
        'N': 'Asparagine',
        'P': 'Proline',
        'Q': 'Glutamine',
        'R': 'Arginine',
        'S': 'Serine',
        'T': 'Threonine',
        'V': 'Valine',
        'W': 'Tryptophan',
        'X': 'Any',
        'Y': 'Tyrosine',
        'Z': '[Glutamine or Glutamic Acid]',
        '*': 'Stop'
    };


    complements.dna = {
        'a': 't',
        'c': 'g',
        'g': 'c',
        't': 'a',
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A',
        'r' : 'y',
        'y' : 'r',
        'R' : 'Y',
        'Y' : 'R',
        'k' : 'm',
        'm' : 'k',
        'K' : 'M',
        'M' : 'K',
        'b' : 'v',
        'v' : 'b',
        'B' : 'V',
        'V' : 'B',
        'd' : 'h',
        'h' : 'd',
        'D' : 'H',
        'H' : 'D',
        'n' : 'n',
        'N' : 'N'
    };
    complements.rna = {
        'a': 'u',
        'c': 'g',
        'g': 'c',
        'u': 'a',
        'A': 'U',
        'C': 'G',
        'G': 'C',
        'U': 'A',
        'r' : 'y',
        'y' : 'r',
        'R' : 'Y',
        'Y' : 'R',
        'k' : 'm',
        'm' : 'k',
        'K' : 'M',
        'M' : 'K',
        'b' : 'v',
        'v' : 'b',
        'B' : 'V',
        'V' : 'B',
        'd' : 'h',
        'h' : 'd',
        'D' : 'H',
        'H' : 'D',
        'n' : 'n',
        'N' : 'N'
    };
    complements.transcribe = {
        'a': 'u',
        'c': 'g',
        'g': 'c',
        't': 'a',
        'A': 'U',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    };
    complements.reverse_transcribe = {
        'a': 't',
        'c': 'g',
        'g': 'c',
        'u': 'a',
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'U': 'A'
    };
    // standard genetic code
    // see https://github.com/cathalgarvey/PySplicer/blob/master/pysplicer/translationtables.py for more
    complements.translate = {
        'uuu' : 'F',    'ucu' : 'S',    'uau' : 'Y',    'ugu' : 'C',
        'uuc' : 'F',    'ucc' : 'S',    'uac' : 'Y',    'ugc' : 'C',
        'uua' : 'L',    'uca' : 'S',    'uaa' : '*',    'uga' : '*',
        'uug' : 'L',    'ucg' : 'S',    'uag' : '*',    'ugg' : 'W',

        'cuu' : 'L',    'ccu' : 'P',    'cau' : 'H',    'cgu' : 'R',
        'cuc' : 'L',    'ccc' : 'P',    'cac' : 'H',    'cgc' : 'R',
        'cua' : 'L',    'cca' : 'P',    'caa' : 'Q',    'cga' : 'R',
        'cug' : 'L',    'ccg' : 'P',    'cag' : 'Q',    'cgg' : 'R',

        'auu' : 'I',    'acu' : 'T',    'aau' : 'N',    'agu' : 'S',
        'auc' : 'I',    'acc' : 'T',    'aac' : 'N',    'agc' : 'S',
        'aua' : 'I',    'aca' : 'T',    'aaa' : 'K',    'aga' : 'R',
        'aug' : 'M',    'acg' : 'T',    'aag' : 'K',    'agg' : 'R',

        'guu' : 'V',    'gcu' : 'A',    'gau' : 'D',    'ggu' : 'G',
        'guc' : 'V',    'gcc' : 'A',    'gac' : 'D',    'ggc' : 'G',
        'gua' : 'V',    'gca' : 'A',    'gaa' : 'E',    'gga' : 'G',
        'gug' : 'V',    'gcg' : 'A',    'gag' : 'E',    'ggg' : 'G'
    };
    complements.reverse_translate = {
        A : ['gcu', 'gcc', 'gca', 'gcg'],
        C : ['ugu', 'ugc'],
        D : ['gau', 'gac'],
        E : ['gaa', 'gag'],
        F : ['uuu', 'uuc'],
        G : ['ggu', 'ggc', 'gga', 'ggg'],
        H : ['cau', 'cac'],
        I : ['auu', 'auc', 'aua'],
        K : ['aaa', 'aag'],
        L : ['uua', 'uug', 'cuu', 'cuc', 'cua', 'cug'],
        M : ['aug'],
        N : ['aau', 'aac'],
        P : ['ccu', 'ccc', 'cca', 'ccg'],
        Q : ['caa', 'cag'],
        R : ['cgu', 'cgc', 'cga', 'cgg', 'aga', 'agg'],
        S : ['ucu', 'ucc', 'uca', 'ucg', 'agu', 'agc'],
        T : ['acu', 'acc', 'aca', 'acg'],
        V : ['guu', 'guc', 'gua', 'gug'],
        W : ['ugg'],
        Y : ['uau', 'uac'],
        '*' : ['uaa', 'uag', 'uga']
    };
    complements.reverse_translate_regexp = {
        A : /gc[acgturyswkmbdhvn]/,
        C : /[tu]g[ctuy]/,
        D : /ga[tcuy]/,
        E : /ga[agr]/,
        F : /[tu][tu][tcuy]/,
        G : /gg[acgturyswkmbdhvn]/,
        H : /ca[tcuy]/,
        I : /a[tu][atcuwmhy]/,
        K : /aa[agr]/,
        L : /c[tu][acgturyswkmbdhvn]|[tu][tu][agr]|[ctuy][tu][agr]/,
        M : /a[tu]g/,
        N : /aa[tucy]/,
        P : /cc[acgturyswkmbdhvn]/,
        Q : /ca[agr]/,
        R : /cg[acgturyswkmbdhvn]|ag[agr]|[cam]g[agr]/,
        S : /[tu]c[acgturyswkmbdhvn]|ag[ct]/,
        T : /ac[acgturyswkmbdhvn]/,
        V : /g[tu][acgturyswkmbdhvn]/,
        W : /[tu]gg/,
        Y : /[tu]a[ctuy]/,
        '*' : /[tu]a[agr]|[tu]ga|[tu][agr]a/
    };


    //todo
    frequencies.ecoli = {

    };


    //deoxy bases e.g. dAMP, dTTP etc. - no salt
    weights.dna = {
        'a' : 313.2,
        't' : 304.2,
        'c' : 289.2,
        'g' : 329.2,
        'u' : 290.2
    };
    //AMP, TMP etc. no salt
    weights.rna = {
        'a' : 329.2,
        't' : 320.2,
        'c' : 305.2,
        'g' : 345.2,
        'u' : 306.2
    };




    var removeNonLetters = function (input) {
        return input.replace(regexps.letters, '');
    };

    var dnaOnly = function (sequence, strict) {
        var filter = strict ? regexps.dna : regexps.dna_degenerate;
        return sequence.replace(filter, '');
    };

    // use regexps.strip.<type> e.g. regexps.strip.dna
    var verify = function(sequence, reg) {
        return (sequence.search(reg) == -1)
    };



    var reverse = function(sequence) {
        return sequence.split("").reverse().join("")
    };

    var complement = function(sequence) {
        return sequence.replace(regexps.dna, function(m) { return complements.dna[m] });
    };

    var rna_to_dna = function (sequence) {
        return sequence.replace(/[uU]/ig, function(m) {return {'u':'t', 'U':'T'}[m]});
    };

    var dna_to_rna = function (sequence) {
        return sequence.replace(/[tT]/ig, function(m) {return {'t':'u', 'T':'U'}[m]});
    };

    var amino_one_to_three = function (sequence) {
        var protein = [];
        angular.forEach(sequence.split(''), function (amino, ind) {protein.push( maps.amino_one_to_three[amino])});
        return protein.join(' ');
    };

    var amino_one_to_full = function (sequence) {
        var protein = [];
        angular.forEach(sequence.split(''), function (amino, ind) {protein.push (maps.amino_one_to_full[amino])});
        return protein.join(' ');
    };

    var revcomp = function(sequence) {
        return reverse(complement(sequence));
    };

    var transcribe = function (sequence) {
        return sequence.replace(regexps.dna, function(m) { return complements.transcribe[m] });
    };

    var reverseTranscribe = function (sequence) {
        return sequence.replace(regexps.rna, function(m) { return complements.reverse_transcribe[m] });
    };

    /**
     * @name DNA.translate
     * @param {string} sequence RNA Sequence to be translated. DNA will be converted (not transcribed) to RNA.
     * @returns {string} Polypeptide translated
     */
    var translate = function (sequence) {
        //todo - multiple orfs instead of just trimming

        if (verify(sequence, regexps.strip.dna)) {
            sequence = dna_to_rna(sequence);
        }

        var seqlen = sequence.length;
        if (seqlen < 3) return '';

        sequence = angular.lowercase(sequence);

        if ((seqlen % 3) != 0)
            sequence = sequence.substring(0, (Math.floor(seqlen / 3) * 3));

        var polypep = "";
        for (var i = 0; i < sequence.length; i = i+3) {
            polypep = polypep + complements.translate[sequence.substr(i, 3)];
        }

        return polypep;
    };

    var reverseTranslate = function (sequence, codonFreq) {
        if (!sequence) return;

        var dumb = (!codonFreq || angular.isUndefined(codonFreq) || codonFreq == null) ? true : false;
        var rna = '';

        for(var seq = sequence.split(''), i = 0; i < seq.length; i++) {
            var options = (complements.reverse_translate[seq[i]]);
            if (dumb) {
                rna = rna + options[Math.floor(Math.random()*options.length)]
            } else {
                //todo implement non-dumb version once have frequency table...
                rna = rna + options[0]
            }
        }

        return rna;
    };

    var randomSequence = function (length, units) {
        units = units || monomers.dna;
        if (typeof length == 'string')
            length = length.length;

        var sequenceArray = [],
            tempNum = 0,
            tempChar = "";
        for (var j = 0; j < (length); j++)	{
            tempNum = Math.floor(Math.random() * units.length);
            tempChar = units[tempNum];
            sequenceArray.push(tempChar);
        }
        return sequenceArray.join("");
    };

    var shuffleSequence = function (sequence) {
        return ($filter('shuffle')(sequence.split(''))).join('');
    };

    /**
     *
     * @param {string} dna Sequence of non-degenerate DNA or RNA
     * @returns {number} fraction of GC in whole sequence
     * @description
     * Doesn't work for degenerate sequences
     */
    var GC_content = function(dna) {
        var full = dna.length,
            gc = (dna.replace(/[^GCS]/gi, '')).length;

        return (gc / full);
    };

    var monomer_count = function (sequence, units) {
        units = (typeof units != 'undefined' ? units : monomers.all);
        var counts = {};

        angular.forEach(units, function(unit) {
            counts[unit] = (sequence.replace(new RegExp('[^' + unit + ']', 'gi'), '')).length;
        });

        return counts;
    };

    var neighbor_count = function (sequence, units) {
        units = (typeof units != 'undefined' ? units : monomers.dna);
        var counts = {};
        angular.forEach(units, function (u1) {
            angular.forEach(units, function (u2) {
                var combo = u1 + u2,
                    matches = sequence.match(new RegExp(combo, 'gi'));
                counts[combo] = matches != null ? matches.length : 0;
            });
        });
        return counts;
    };

    //todo calc: rlnk, deltaG, H, S @ standard cond
    var gibbs = function(sequence) {
        var dg = 0,
            dh = 0,
            ds = 0,
            counts = monomer_count(sequence, monomers.nucleotide),
            neighbors = neighbor_count(sequence, monomers.nucleotide),
            gc = GC_content(sequence);

        // convert to moles and then adjust for greater stabilizing effects of Mg compared to Na or K.
        //salt+=(mg/Math.pow(10,3) * 140);

        return [dg, dh, ds]

    };



    //todo - more rigorous - account for primer concencration, pH, mg -- should rely on Gibbs
    //todo - write for RNA (need to account for folding etc.)
    var melting_temp = function(dna, conc) {
        if (dna.length < 1) return;

        if ((regexps.strip.dna).exec(dna)) {
            console.log('contains non-dna letters... aborting melting temp calc');
            return;
        }

        conc = typeof conc != 'undefined' ? conc : {};
        conc.dna = typeof conc.dna != 'undefined' ? conc.dna : 0.050;     // Molar
        conc.salt = typeof conc.salt != 'undefined' ? conc.salt : 0.050;   // Molar
        conc.mg = typeof conc.mg != 'undefined' ? conc.mg : 0.0025;        // Molar
        conc.ph = typeof conc.ph != 'undefined' ? conc.ph : 7.0;        //pH

        var counts = monomer_count(dna, monomers.dna);

        if (dna.length < 14) {
            return (counts['a']+counts['t'])*2 + (counts['c'] + counts['g'])*4 - 16.6*(Math.log(0.050) / Math.log(10)) + (16.6 * (Math.log(conc.salt) / Math.log(10)))
        } else if (dna.length < 50) {
            return 100.5 + (41 * (counts['g']+counts['c'])/(counts['a']+counts['t']+counts['g']+counts['c'])) - (820/(counts['a']+counts['t']+counts['g']+counts['c'])) + (16.6 * (Math.log(conc.salt) / Math.log(10)))
        } else {
            return 81.5 + (41 * (counts['g']+counts['c'])/(counts['a']+counts['t']+counts['g']+counts['c'])) - (500/(counts['a']+counts['t']+counts['g']+counts['c'])) + (16.6 * (Math.log(conc.salt) / Math.log(10)))
        }

    };

    //todo - handle degenerate, handle proteins, verify correct
    var molecular_weight = function(sequence) {
        if (!sequence || sequence.length < 1) return;

        var type = verify(sequence, regexps.strip.dna) ? 'dna' : 'rna',
            nucs = monomers.nucleotide,
            weight = 0,
            counts = monomer_count(sequence, nucs);

        angular.forEach(nucs, function(nuc) {
            weight = weight + (weights[type][nuc] * counts[nuc]);
        });

        // future - more rigorous
        //triphospate --- could substract 61.96 instead to remove phospate and add hydroxyl, and subsequently can add 79.0 for restriction-enzyme cut DNA which leaves 5' monophospate
        weight = weight + 159.0;

        return weight;
    };


    return {
        regexps: regexps,
        monomers : monomers,
        maps : maps,
        complements : complements,
        frequencies : frequencies,
        weights : weights,

        removeNonLetters : removeNonLetters,
        dnaOnly : dnaOnly,
        verify : verify,

        reverse : reverse,
        complement : complement,
        rna_to_dna : rna_to_dna,
        dna_to_rna : dna_to_rna,
        amino_one_to_three : amino_one_to_three,
        amino_one_to_full : amino_one_to_full,
        revcomp : revcomp,
        transcribe : transcribe,
        reverseTranscribe : reverseTranscribe,
        translate : translate,
        reverseTranslate : reverseTranslate,

        randomSequence : randomSequence,
        shuffleSequence : shuffleSequence,
        GC_content : GC_content,
        monomer_count : monomer_count,
        neighbor_count : neighbor_count,
        gibbs : gibbs,
        melting_temp : melting_temp,
        molecular_weight : molecular_weight
    }
}]);

Application.Dna.controller('dnaCtrl', ['$scope', 'DNA', function($scope, DNA) {
    $scope.DNA = DNA;

    //todo - if possible, assign these directly from DOM using a directive...
    $scope.$watch('sequence', function(newval, oldval) {
        $scope.rna = DNA.transcribe(newval);
    });

    $scope.$watch('rna', function (newval, oldval) {
        $scope.protein = DNA.translate(newval);
    });

    $scope.sequence = 'aaatttgggcccatgcta';
}]);

Application.Dna.directive('dnaShuffle', ['DNA', function(DNA) {
    return {
        restrict: 'A',
        require: 'ngModel',
        link: function(scope, element, attr, ngModel) {
            ngModel.$formatters.push(DNA.shuffleSequence)
        }
    };
}]);

Application.Dna.directive('dnaRandom', ['DNA', function(DNA) {
    return {
        restrict: 'A',
        require: 'ngModel',
        link: function(scope, element, attr, ngModel) {
            ngModel.$formatters.push(DNA.randomSequence)
        }
    };
}]);

Application.Dna.directive('dnaReverse', ['DNA', function(DNA) {
    return {
        restrict: 'A',
        require: 'ngModel',
        link: function(scope, element, attr, ngModel) {
            ngModel.$formatters.push(DNA.reverse)
        }
    };
}]);

Application.Dna.directive('dnaComplement', ['DNA', function(DNA) {
    return {
        restrict: 'A',
        require: 'ngModel',
        link: function(scope, element, attr, ngModel) {
            ngModel.$formatters.push(DNA.complement)
        }
    };
}]);

Application.Dna.directive('dnaRevcomp', ['DNA', function(DNA) {
    return {
        restrict: 'A',
        require: 'ngModel',
        link: function(scope, element, attr, ngModel) {
            ngModel.$formatters.push(DNA.revcomp)
        }
    };
}]);

Application.Dna.directive('dnaTranscribe', ['DNA', function(DNA) {
    return {
        restrict: 'A',
        require: 'ngModel',
        link: function(scope, element, attr, ngModel) {

            var fn = function(input) {
                console.log(input);
                return DNA.transcribe(input);
            };

            ngModel.$formatters.push(fn)
        }
    };
}]);

Application.Dna.directive('dnaReverseTranscribe', ['DNA', function(DNA) {
    return {
        restrict: 'A',
        require: 'ngModel',
        link: function(scope, element, attr, ngModel) {

            var fn = function(input) {
                return DNA.reverseTranscribe(input);
            };

            ngModel.$formatters.push(fn)
        }
    };
}]);

Application.Dna.directive('dnaTranslate', ['DNA', function(DNA) {
    return {
        restrict: 'A',
        require: 'ngModel',
        link: function(scope, element, attr, ngModel) {
            var fn = function(input) {
                return DNA.translate(input);
            };

            ngModel.$formatters.push( fn );
        }
    };
}]);

Application.Dna.directive('dnaReverseTranslate', ['DNA', '$parse', function(DNA, $parse) {
    return {
        restrict: 'A',
        require: 'ngModel',
        link: function(scope, element, attr, ngModel) {
            //if want to pass options to function
            //var flags = JSON.parse(attr.dnaReverseTranslate);
            //console.log(flags);

            var codonFreq = $parse(attr.dnaReverseTranslate)(scope) || false;

            var fn = function(input) {
                return DNA.reverseTranslate(input, codonFreq);
            };

            ngModel.$formatters.push( fn );
        }
    };
}]);

Application.Dna.directive('dnaDnaOnly', ['DNA', function(DNA) {
    return {
        restrict: 'A',
        require: 'ngModel',
        link: function(scope, element, attr, ngModel) {
            ngModel.$formatters.push(DNA.dnaOnly)
        }
    };
}]);

Application.Dna.directive('dnaRemoveNonLetters', ['DNA', function(DNA) {
    return {
        restrict: 'A',
        require: 'ngModel',
        link: function(scope, element, attr, ngModel) {
            ngModel.$formatters.push(DNA.removeNonLetters)
        }
    };
}]);

Application.Dna.directive('dnaAminoOneToThree', ['DNA', function(DNA) {
    return {
        restrict: 'A',
        require: 'ngModel',
        link: function(scope, element, attr, ngModel) {
            ngModel.$formatters.push(DNA.amino_one_to_three)
        }
    };
}]);

Application.Dna.directive('dnaAminoOneToFull', ['DNA', function(DNA) {
    return {
        restrict: 'A',
        require: 'ngModel',
        link: function(scope, element, attr, ngModel) {
            ngModel.$formatters.push(DNA.amino_one_to_full)
        }
    };
}]);

Application.Dna.directive('dnaGcContent', ['DNA', '$filter', function(DNA, $filter) {
    return {
        restrict: 'A',
        require: 'ngModel',
        link: function(scope, element, attr, ngModel) {
            var fn = function (input) {
                return $filter('number')(DNA.GC_content(input))
            };
            ngModel.$formatters.push(fn)
        }
    };
}]);

Application.Dna.directive('dnaMeltingTemp', ['DNA', '$filter', function(DNA, $filter) {
    return {
        restrict: 'A',
        require: 'ngModel',
        link: function(scope, element, attr, ngModel) {
            var fn = function (input) {
                return $filter('number')(DNA.melting_temp(input))
            };
            ngModel.$formatters.push(fn)
        }
    };
}]);


Application.Dna.directive('dnaMolecularWeight', ['DNA', '$filter', function(DNA, $filter) {
    return {
        restrict: 'A',
        require: 'ngModel',
        link: function(scope, element, attr, ngModel) {
            var fn = function (input) {
                return $filter('number')(DNA.molecular_weight(input))
            };
            ngModel.$formatters.push(fn)
        }
    };
}]);

//testing only
Application.Dna.directive('dnaGibbs', ['DNA', function(DNA) {
    return {
        restrict: 'A',
        require: 'ngModel',
        link: function(scope, element, attr, ngModel) {
            var fn = function (input) {
                return (DNA.gibbs(input))[0];
            };
            ngModel.$formatters.push(fn)
        }
    };
}]);