'use strict';

//todo - functionality
// directives that transform don't affect ngModel on change
// codon frequency e.g. e-coli --->> smart backTranslate
// gibbs calculation
// restriction site service
// FASTA reader, genbank
// alignment
// mutation
// pcr products - annealing etc.
// fuzzy search

Application.Dna.service('DNA', [function() {

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
        geneticCodes = {},
        frequencies = {},
        entropies = {},
        weights = {};

    regexps.letters = /[A-Z]/ig;
    regexps.nucleotide = /[gautcGAUTC]/g;
    regexps.dna = /[gatcGATC]/g;
    regexps.rna = /[gaucGAUC]/g;
    regexps.nucleotide_degnerate = /[gatucryswkmbdhvnxGATUCRYSWKMBDHVNX\^_]/g;
    regexps.protein = /[ACDEFGHIKLMNPQRSTVWYZacdefghiklmnpqrstvwyz\*]/g;
    regexps.protein_degenerate = /[ABCDEFGHIKLMNPQRSTVWYXZabcdefghiklmnpqrstvwyxz\*]/g;

    regexps.strip = {};
    regexps.strip.dna = /[^gatcGATC]/g;
    regexps.strip.nucleotide_degenerate = /[^gatucryswkmbdhvnxGATUCRYSWKMBDHVNX]/g;
    regexps.strip.rna = /[^gaucGAUC]/g;
    regexps.strip.protein = /[^ACDEFGHIKLMNPQRSTVWYZacdefghiklmnpqrstvwyz\*]/g;
    regexps.strip.protein_degenerate = /[^ABCDEFGHIKLMNPQRSTVWYXZabcdefghiklmnpqrstvwyxz\*]/g;
    regexps.strip.alignment = /[^\.\-]/g;
    regexps.strip.non_letters = /[^A-Z]/ig;
    regexps.strip.whitespace = /\s/g;


    monomers.dna = ('acgt').split('');
    monomers.rna = ('acug').split('');
    monomers.nucleotide = ('acgtu').split('');
    monomers.nucleotide_degenerate = ('acgturyswkmbdhvnx').split('');
    monomers.protein = ('ACDEFGHIKLMNPQRSTVWYZ*').split('');
    monomers.protein_degnerate = ('ABCDEFGHIKLMNPQRSTVWXYZ*').split('');
    //dna, rna, protein
    monomers.all = ("ABCDEFGHIKLMNPQRSTUVWXYZ*").split('');


    maps.nucleotide_degenerate = {
        'A': 'A', 'B': '[CGTU]', 'C': 'C', 'D': '[AGTU]', 'G': 'G', 'H': '[ACTU]',
        'K': '[GTU]', 'M': '[AC]', 'N': '[ACGTU]', 'R': '[AG]', 'S': '[CG]',
        'T': 'T', 'U': 'U', 'V': '[ACG]', 'W': '[ATU]', 'Y': '[CTU]',
        '.': '[ACGTU]', '-': '[ACGTU]'};
    maps.nucleotide_undegenerate = {A: "A", C: "C", G: "G", T: "T", U: "U", '[ACGTU]': "N", '[ACG]': "V", '[ACTU]': "H", '[AC]': "M", '[AGTU]': "D", '[AG]': "R", '[ATU]': "W", '[CGTU]': "B", '[CG]': "S", '[CTU]': "Y", '[GTU]': "K"};
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
        'N' : 'N',
        '^' : '_',
        '_' : '^'
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
        'N' : 'N',
        '^' : '_',
        '_' : '^'
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
    //https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c
    geneticCodes.standard = {};



    //todo
    //see
    frequencies.ecoli = {

    };


    //used in melting point calculations -- currently only DNA
    //SantaLucia, J. (1998) Proc. Nat. Acad. Sci. USA 95, 1460.
    entropies.ds_terminal = {
        g : -2.8,
        a : 4.1,
        t : 4.1,
        c : -2.8
    };
    entropies.dh_terminal = {
        g : 0.1,
        a : 2.3,
        t : 2.3,
        c : 0.1
    };
    entropies.ds = {
        gg: -19.9, ag: -21.0, tg: -22.7, cg: -27.2,
        ga: -22.2, aa: -22.2, ta: -21.3, ca: -22.7,
        gt: -22.4, at: -20.4, tt: -22.2, ct: -21.0,
        gc: -27.2, ac: -22.4, tc: -22.2, cc: -19.9
    };
    entropies.dh = {
        gg : -8.0, ag : -7.8, tg : -8.5, cg : -10.6,
        ga : -8.2, aa : -7.9, ta : -7.2, ca : -8.5,
        gt : -8.4, at : -7.2, tt : -7.9, ct : -7.8,
        gc : -10.6,ac : -8.4, tc : -8.2, cc : -8.0
    };


    //verify
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



    /**************
    Basics
    **************/

    var lettersOnly = function (input) {
        return input.replace(regexps.letters, '');
    };

    var dnaOnly = function (sequence, strict) {
        var filter = strict ? regexps.strip.dna : regexps.strip.nucleotide_degenerate;
        return sequence.replace(filter, '');
    };

    // use regexps.strip.<type> e.g. regexps.strip.dna
    var verify = function(sequence, reg) {
        return (sequence.search(reg) == -1)
    };

    //given a sequence (e.g. ACNT), return RegExp friendly version using given map (e.g. AC[ACGTU]T)
    //note currently only nucleotides
    var undegenerize = function(sequence, map) {
        map = !!map ? map : maps.nucleotide_degenerate;
        return sequence.replace(regexps.nucleotide_degenerate, function(m) {return map[m.toUpperCase()]});
    };


    /**************
    Sequence manipulation
    **************/

    var reverse = function(sequence) {
        return sequence.split("").reverse().join("")
    };

    var complement = function(sequence) {
        return sequence.replace(regexps.nucleotide_degnerate, function(m) {return complements.dna[m] } );
    };

    var rna_to_dna = function (sequence) {
        return sequence.replace(/[uU]/ig, function(m) {return {'u':'t', 'U':'T'}[m]});
    };

    var dna_to_rna = function (sequence) {
        return sequence.replace(/[tT]/ig, function(m) {return {'t':'u', 'T':'U'}[m]});
    };

    var amino_one_to_three = function (sequence) {
        var protein = [];
        _.forEach(sequence.split(''), function (amino, ind) {protein.push( maps.amino_one_to_three[amino])});
        return protein.join(' ');
    };

    var amino_one_to_full = function (sequence) {
        var protein = [];
        _.forEach(sequence.split(''), function (amino, ind) {protein.push (maps.amino_one_to_full[amino])});
        return protein.join(' ');
    };

    var revcomp = function(sequence) {
        return reverse(complement(sequence));
    };

    var randomSequence = function (length, units) {
        units = units || monomers.dna;
        if (typeof length == 'string')
            length = length.length;

        if (length == 0) return '';

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
        return _.shuffle(sequence.split('')).join('');
    };

    //adds a space every three characters
    var codonSpace = function (sequence) {
        return sequence.match(/.{1,3}/g).join(" ");
    };

    /**************
     Central Dogma
     **************/

    var transcribe = function (sequence) {
        return sequence.replace(regexps.dna, function(m) { return complements.transcribe[m] });
    };

    var reverseTranscribe = function (sequence) {
        return sequence.replace(regexps.rna, function(m) { return complements.reverse_transcribe[m] });
    };

    /**
     * @name DNA.translate
     * @param {string} sequence RNA Sequence to be translated. DNA will be converted (not transcribed) to RNA.
     * @param {number} forceOffset Force an offset of 0,1,2 bases as reading frame
     * @returns {string} Polypeptide translated
     */
    var translate = function (sequence, forceOffset) {

        //checks
        if (verify(sequence, regexps.strip.dna)) {
            sequence = dna_to_rna(sequence);
        }
        var seqlen = sequence.length;
        if (seqlen < 3) return '';
        sequence = sequence.toLowerCase();


        var offset = (!!forceOffset) ? forceOffset : probableORF(sequence),
            polypep = "";

        sequence = sequence.substring(offset);

        if ((seqlen % 3) != 0)
            sequence = sequence.substring(0, (Math.floor(seqlen / 3) * 3));

        for (var i = 0; i < sequence.length; i = i+3) {
            polypep = polypep + complements.translate[sequence.substr(i, 3)];
        }

        return polypep;
    };

    var reverseTranslate = function (sequence, codonFreq) {
        if (!sequence) return;

        var dumb = (!codonFreq || _.isUndefined(codonFreq) || codonFreq == null) ? true : false;
        var rna = '';

        for(var seq = sequence.split(''), i = 0; i < seq.length; i++) {
            var options = (complements.reverse_translate[seq[i]]);
            if (dumb) {
                rna = rna + options[Math.floor(Math.random()*options.length)]
            } else {
                //future implement non-dumb version once have frequency table...
                rna = rna + options[0]
            }
        }

        return rna;
    };


    /**
     * @description Determines fragments for each ORF for a given sequence
     * @param sequence
     * @returns {boolean|object} false if no proteins (atg...stop), otherwise object in the following form, with frags descending by length:
     {
        sequence : <input sequence>,
        frags : [
            {
                match : <matched seq (from atg ... tag|tga|taa)>
                index : <match index>
                length : <dna seq length>
            }
            ...
        ],
        longest : <length of longest putative protein>
     }

     */
    var calcORFs = function(sequence) {
        var frames = {},
            orfReg = new RegExp('((atg)(.+?)(ta[agr]|tga|t[agr]a))', 'gi');

        //no real proteins
        if (!orfReg.test(sequence)) return false;

        for (var i = 0; i < 3; i++) {
            frames[i] = {};
            frames[i].sequence = sequence;
            frames[i].frags = [];

            var match;
            while ((match = orfReg.exec(sequence)) != null) {
                frames[i].frags.push({
                    match : match[0],
                    index : match.index,
                    length : match[0].length
                })
            }
            frames[i].frags.sort(function(a, b){ return b.length - a.length; });
            frames[i].longest = frames[i].frags[0].length;
        }

        return frames;
    };

    /**
     * @description Returns the offset for the longest ORF
     * @param sequence
     * @returns {number}
     */
    var probableORF = function (sequence) {
        var offset = 0,
            longest = 0;

        var frames = calcORFs(sequence);

        //no proteins, retruned false
        if (!frames) return 0;

        for (var i = 0; i < frames.length; i++) {
            if (longest < frames[i].longest) {
                longest = frames[i].longest;
                offset = i;
            }
        }
        //future - check codon frequency
        //future - alignment to database e.g. blast or something

        return offset;
    };


    /**************
     Quantification
     **************/

    /**
     * @description Returns number of occurances of a given base or oligo in a sequence. Handles overlaps.
     * @example occuranceCount('aataa', 'a') -> 4
     * @example occuranceCount('aaaa', 'aa') -> 3
     * @param sequence
     * @param {String} oligo Single-letter base or oligo to look for
     * @param {boolean=} convertDegenerate Whether degenerate nucleotides should be converted to standard [ACGTU]. Default is false
     * @returns {Number} Number of matches found, 0 if none found
     */
    var occuranceCount = function(sequence, oligo, convertDegenerate) {
        var search = !!convertDegenerate ? undegenerize(oligo) : oligo;
        var matches = sequence.match(new RegExp('(?=(' + search + '))', 'gi'));
        return !!matches ? matches.length : 0;
    };

    /**
     * @description Get a count of each monomer in a sequence
     * @param {string} sequence
     * @param {object=} units Array or Object of monomers to search e.g. monomers.nucleotides
     * @returns {object} object whose keys are monomers, and values are number of occurances
     */
    var monomer_count = function (sequence, units) {
        units = (typeof units != 'undefined' ? units : monomers.all);
        var counts = {};
        _.forEach(units, function(unit) {
            counts[unit] = occuranceCount(sequence, unit);
        });
        return counts;
    };

    /**
     * @description Determines number of occurances of neighboring nucleotides, e.g. 'ac' -> 4
     * @note Currently only DNA
     * @param sequence
     * @param minCount
     * @param units
     * @returns {object}
     */
    var neighbor_count = function (sequence, minCount, units) {
        units = (typeof units != 'undefined') ? units : monomers.dna;
        minCount = (typeof minCount != 'undefined') ? minCount : 0;
        var counts = {};
        _.forEach(units, function (u1) {
            _.forEach(units, function (u2) {
                var combo = u1 + u2,
                    count = occuranceCount(sequence, combo);
                if (count >= minCount)
                    counts[combo] = count;
            });
        });
        return counts;
    };

    /**
     * @description Calculates GC content of non-degenerate nucleotide sequence
     * @note Doesn't work for degenerate sequences
     * @param {string} dna Sequence of non-degenerate DNA or RNA
     * @returns {number} fraction of GC in whole sequence
     */
    var GC_content = function(dna) {
        var full = dna.length,
            gc = (dna.replace(/[^GCS]/gi, '')).length;

        return (gc / full);
    };


    var gibbs = function(sequence) {
        var dg = 0,
            dh = 0,
            ds = 0,
            counts = monomer_count(sequence, monomers.nucleotide),
            neighbors = neighbor_count(sequence, monomers.nucleotide),
            gc = GC_content(sequence);

        //todo calc: rlnk, deltaG, H, S @ standard cond

        return [dg, dh, ds]

    };


    var melting_temp_basic = function (sequence) {
        if (sequence.length < 14) {
            var counts = monomer_count(sequence, monomers.dna);
            return (counts['a']+counts['t'])*2 + (counts['c'] + counts['g'])*4
        }

        //64.9°C + 41°C x (number of G’s and C’s in the oligo – 16.4)/N
        return 64.9+41*(GC_content(sequence)*sequence.length - 16.4)/sequence.length;
    };

    //verify
    var melting_temp_saltAdjusted = function (sequence, saltConc) {
        saltConc = !!saltConc ? saltConc : 0.050;

        var counts = monomer_count(sequence, monomers.dna);

        if (sequence.length < 14) {
            return (counts['a']+counts['t'])*2 + (counts['c'] + counts['g'])*4 - 16.6*(Math.log(0.050) / Math.log(10)) + (16.6 * (Math.log(saltConc) / Math.log(10)))
        } else if (sequence.length < 50) {
            return 100.5 + (41 * (counts['g']+counts['c'])/(counts['a']+counts['t']+counts['g']+counts['c'])) - (820/(counts['a']+counts['t']+counts['g']+counts['c'])) + (16.6 * (Math.log(saltConc) / Math.log(10)))
        } else {
            return 81.5 + (41 * (counts['g']+counts['c'])/(counts['a']+counts['t']+counts['g']+counts['c'])) - (500/(counts['a']+counts['t']+counts['g']+counts['c'])) + (16.6 * (Math.log(saltConc) / Math.log(10)))
        }
    };

    /**
     * @description More complex melting temperature calculation for a given sequence. The most sophisticated Tm calculations take into account the exact sequence and base stacking parameters, not just the base composition.
     * Tm = ((1000* dh)/(ds+(R * Math.log(primer concentration))))-273.15;
     * @param {string} sequence ONLY CANONICAL DNA BASES (A,C,G,T)
     * @param {object} conc Concentrations with parameters 'dna' 'salt' and 'mg' all Molar
     */

    //todo - account for buffers, see NEB site
    //https://www.neb.com/tools-and-resources/interactive-tools/tm-calculator
    var melting_temp = function (sequence, conc) {
        if (sequence.length < 1) return;

        if ((regexps.strip.dna).exec(sequence)) {
            console.log('contains non-dna letters... aborting melting temp calc');
            return;
        }

        sequence = sequence.toLowerCase();

        conc = typeof conc != 'undefined' ? conc : {};
        conc.dna = typeof conc.dna != 'undefined' ? conc.dna : 0.0000002;   // Molar (200nM)
        conc.salt = typeof conc.salt != 'undefined' ? conc.salt : 0.050;    // Molar (50mM)
        conc.mg = typeof conc.mg != 'undefined' ? conc.mg : 0.0015;         // Molar (2.5mM)
        //conc.ph = typeof conc.ph != 'undefined' ? conc.ph : 7.0;          //pH - NOT USED CURRENTLY

        var R = 1.987; //universal gas constant in Cal/degrees C * mol
        var ds = 0;    //cal/Kelvin/mol
        var dh = 0;    //kcal/mol

        // salt correction
        var correctedSalt = conc.salt + conc.mg * 140; //adjust for greater stabilizing effects of Mg compared to Na or K. See von Ahsen et al 1999
        ds = ds + 0.368 * (sequence.length - 1) * Math.log(correctedSalt); //from von Ahsen et al 1999

        // terminal corrections
        ds = ds + entropies.ds_terminal[sequence.charAt(0)] + entropies.ds_terminal[sequence.charAt(sequence.length - 1)];
        dh = dh + entropies.dh_terminal[sequence.charAt(0)] + entropies.dh_terminal[sequence.charAt(sequence.length - 1)];

        //nearest neighbors
        var neighbors = neighbor_count(sequence);
        _.forEach(neighbors, function (num, pair) {
            ds = ds + (entropies.ds[pair] * num);
            dh = dh + (entropies.dh[pair] * num);
        });

        return ((1000 * dh) / (ds + (R * Math.log(conc.dna / 2)))) - 273.15;

    };

    //todo - handle proteins, verify correct
    var molecular_weight = function(sequence) {
        if (!sequence || sequence.length < 1) return;

        var type = verify(sequence, regexps.strip.dna) ? 'dna' : 'rna',
            nucs = monomers.nucleotide,
            weight = 0,
            counts = monomer_count(sequence, nucs);

        _.forEach(nucs, function(nuc) {
            weight = weight + (weights[type][nuc] * counts[nuc]);
        });

        // future - more rigorous
        //triphospate --- could substract 61.96 instead to remove phospate and add hydroxyl, and subsequently can add 79.0 for restriction-enzyme cut DNA which leaves 5' monophospate
        weight = weight + 159.0;

        return weight;
    };


    var selfComplimentarity = function (rna) {
        //todo calculate 2° structure
    };



    /**************
     Species Specific
     **************/

    //todo
    var codonOptimize = function(sequence, species) {

    };

    //todo
    var findSilentSites = function(sequence, species) {

    };



    /**************
     Checks
     **************/

    /**
     * @description
     * @param unit
     * @param repeat
     * @returns {string}
     */
    var createRun = function(unit, repeat) {
        return new Array(+repeat + 1).join(unit);
    };

    /**
     * @description Look for double repeats like 'acacacac'
     * @note currently DNA only
     * @param oligo
     * @param max
     */
    var verifyRepeats = function(oligo, max) {
        var flag = true;
        max = typeof max != 'undefined' ? max : 4;
        monomers = typeof monomers != 'undefined' ? monomers : monomers.nucleotide;

        //first, neighbor count for possible problems
        var counts = neighbor_count(oligo, max);

        //then go through neighbors
        for (var index in counts) {
            if (!counts.hasOwnProperty(index))
                break;
            if ((new RegExp(createRun(index, max), 'ig')).test(oligo))
                flag = false;
        }

        return flag;
    };

    /**
     * @description checks for runs of the same nucleotide (e.g. 'aaaa')
     * @param oligo
     * @param max
     * @param monomers
     * @returns {boolean}
     */
    var verifyRuns = function(oligo, max) {
        var flag = true;
        max = typeof max != 'undefined' ? max : 4;
        monomers = monomers.nucleotide;

        for (var i = 0; i < monomers.length; i++) {
            if ( occuranceCount( oligo, createRun(monomers[i], max) ) )
                flag = false;
        }

        return flag;
    };

    var verifyMeltingTemps = function(primers) {
        var max = 0, min = 0;
        for (var i = 0; i < primers.length; i++) {
            var temp = melting_temp(primers[i]);
            max = (max > temp) ? max : temp;
            min = (min < temp) ? min : temp;
        }

        if (max == 0 || min == 0) {}

        return ((max - min) <= 6);
    };



    /****************
     PARSING -- future
     ****************/

    var parseFASTA = {};

    var parseEMBL = {};

    var parseGenbank = {};

    var outputFASTA = {};

    var outputEMBL = {};

    var outputGenbank = {};


    return {
        //config, info
        regexps: regexps,
        monomers : monomers,
        maps : maps,
        complements : complements,
        geneticCodes : geneticCodes,
        frequencies : frequencies,
        weights : weights,

        //utility
        lettersOnly : lettersOnly,
        dnaOnly : dnaOnly,
        verify : verify,
        undegenerize : undegenerize,

        //manipulation
        reverse : reverse,
        complement : complement,
        rna_to_dna : rna_to_dna,
        dna_to_rna : dna_to_rna,
        amino_one_to_three : amino_one_to_three,
        amino_one_to_full : amino_one_to_full,
        revcomp : revcomp,
        randomSequence : randomSequence,
        shuffleSequence : shuffleSequence,
        codonSpace : codonSpace,

        //central dogma
        transcribe : transcribe,
        reverseTranscribe : reverseTranscribe,
        translate : translate,
        reverseTranslate : reverseTranslate,
        calcORFs : calcORFs,
        probableORF : probableORF,

        //quantification
        occuranceCount : occuranceCount,
        monomer_count : monomer_count,
        neighbor_count : neighbor_count,
        GC_content : GC_content,
        gibbs : gibbs,
        melting_temp_basic : melting_temp_basic,
        melting_temp_saltAdjusted : melting_temp_saltAdjusted,
        melting_temp : melting_temp,
        molecular_weight : molecular_weight,

        //cehcks
        createRun : createRun,
        verifyRepeats : verifyRepeats,
        verifyRuns : verifyRuns,
        verifyMeltingTemps : verifyMeltingTemps
    }
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
            ngModel.$formatters.push(DNA.lettersOnly)
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