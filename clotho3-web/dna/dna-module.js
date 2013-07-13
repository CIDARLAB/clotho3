'use strict';

//todo -
// restriction site service
// FASTA reader
// 1 - to - 3 protein code
// alignment
// mutation
// pcr products
// fuzzy search

Application.Dna.controller('dnaCtrl', ['$scope', 'DNA', function($scope, DNA) {
    $scope.DNA = DNA;

    $scope.sequence = 'actgatctgactgt';

    $scope.complement = DNA.complement($scope.sequence);
    $scope.revcomp = DNA.revcomp($scope.sequence);
}]);

Application.Dna.service('DNA', [function() {

/*
    iupac_nucleotides = [ 'A', 'T', 'C', 'G', 'U',    Canonical bases
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
    frequencies = {};

    regexps.letters = /[A-Z]/ig;
    regexps.nucleotide = /[gautcGAUTC]/g;
    regexps.dna = /[gatcGATC]/g;
    regexps.rna = /[gaucGAUC]/g;
    regexps.dna_degnerate = /[gatucryswkmbdhvnxGATUCRYSWKMBDHVNX]/g;
    regexps.protein = /[ACDEFGHIKLMNPQRSTVWYZacdefghiklmnpqrstvwyz\*]/g;
    regexps.protein_degenerate = /[ABCDEFGHIKLMNPQRSTVWYXZabcdefghiklmnpqrstvwyxz\*]/g;

    regexps.strip = {};
    regexps.strip.dna_strict = /[^gatcGATC]/g;
    regexps.strip.dna = /[^gatucryswkmbdhvnxGATUCRYSWKMBDHVNX]/g;
    regexps.strip.rna = /[^gaucGAUC]/g;
    regexps.strip.protein = /[^ACDEFGHIKLMNPQRSTVWYZacdefghiklmnpqrstvwyz\*]/g;
    regexps.strip.protein_degenerate = /[^ABCDEFGHIKLMNPQRSTVWYXZabcdefghiklmnpqrstvwyxz\*]/g;
    regexps.strip.alignment = /[^\.\-]/g;
    regexps.strip.non_letters = /[^A-Z]/ig;
    regexps.strip.whitespace = /\s/g;


    monomers.dna = ['a', 'c', 't', 'g'];
    monomers.rna = ['a', 'c', 'u', 'g'];
    monomers.dna_degenerate = ['g','a','t','u','c','r','y','s','w','k','m','b','d','h','v','n','x'];
    monomers.protein = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','Z','*'];
    monomers.protein_degnerate = ['A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z','*'];


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
        'B': 'Aspartic Acid or Asparagine',
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
        'Z': 'Glutamine or Glutamic Acid',
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

    var revcomp = function(sequence) {
        return reverse(complement(sequence));
    };

    var transcribe = function (sequence, reverse) {
        var dict = reverse ? complements.reverse_transcribe : complements.transcribe;
        return sequence.replace(regexps.dna, function(m) { return dict[m] });
    };

    var translate = function (sequence, reverse, dumb) {
        //todo - multiple orfs instead of just trimming



        if (verify(sequence, regexps.strip.dna)) {
            sequence = dna_to_rna(sequence);
        }

        var seqlen = sequence.length;
        if (seqlen < 3) return;

        sequence = angular.lowercase(sequence);

        if ((seqlen % 3) != 0) {
            console.log( (Math.floor(seqlen / 3) * 3) -1 );
            sequence = sequence.substring(0, (Math.floor(seqlen / 3) * 3) )
        }

        console.log(sequence);

        var polypep = "";
        for (var i = 0; i < sequence.length; i = i+3) {
            console.log(sequence.substr(i, 3));
            polypep = polypep + complements.translate[sequence.substr(i, 3)];
        }
         
        console.log(polypep);

        return polypep;
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

    var randomSequence = function (length, monomers) {
        monomers = monomers || monomers.dna;
        var sequenceArray = [],
            tempNum = 0,
            tempChar = "";
        for (var j = 0; j < (length); j++)	{
            tempNum = Math.floor(Math.random() * monomers.length);
            tempChar = monomers[tempNum];
            sequenceArray.push(tempChar);
        }
        return sequenceArray.join("");
    };

    var shuffleSequence = function (sequence) {
        var o = sequence.split('');
        for(var j, x, i = o.length; i; j = parseInt(Math.random() * i), x = o[--i], o[i] = o[j], o[j] = x);
        o.join('');
        return o;
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
            gc = (dna.replace(/[atuATU]/, '')).length;

        return (gc / full);
    };


    var melting_temp = function(dna, conc) {
        conc.dna = conc.dna || 400;     // nM
        conc.salt = conc.salt || 50;    // mM
        conc.mg = conc.mg || 2.5;       // mM

        //todo

    };

    var gibbs = function(sequence) {
        //todo calc: rlnk, deltaG, H, S @ standard cond
    };




    return {
        regexps: regexps,
        monomers : monomers,
        maps : maps,
        complements : complements,

        reverse : reverse,
        complement : complement,
        rna_to_dna : rna_to_dna,
        dna_to_rna : dna_to_rna,
        revcomp : revcomp,
        transcribe : transcribe,
        translate : translate,
        removeNonLetters : removeNonLetters,
        dnaOnly : dnaOnly,
        verify : verify,
        randomSequence : randomSequence,
        shuffleSequence : shuffleSequence,
        GC_content : GC_content,
        melting_temp : melting_temp,
        gibbs : gibbs
    }
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

//todo - add ability to reverse, etc.
Application.Dna.directive('dnaTranscribe', ['DNA', function(DNA) {
    return {
        restrict: 'A',
        require: 'ngModel',
        link: function(scope, element, attr, ngModel) {
            ngModel.$formatters.push(DNA.transcribe)
        }
    };
}]);

//todo - add ability to reverse, etc.
Application.Dna.directive('dnaTranslate', ['DNA', function(DNA) {
    return {
        restrict: 'A',
        require: 'ngModel',
        link: function(scope, element, attr, ngModel) {
            ngModel.$formatters.push(DNA.translate)
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
