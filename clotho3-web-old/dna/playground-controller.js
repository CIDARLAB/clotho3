'use strict';

Application.Dna.controller('constructionCtrl', ['$scope', '$parse', 'Clotho', 'DNA', 'Digest', 'PCR', 'Construction', '$http', '$compile', function($scope, $parse, Clotho, DNA, Digest, PCR, Construction, $http, $compile) {
    $scope.DNA = DNA;
    $scope.PCR = PCR;
    $scope.Digest = Digest;
    $scope.Construction = Construction;


    //testing
    var seq = '^acct_acactacgactatcgtagcatgc^ttgg_';
    console.log(DNA.complement(seq));
    console.log(DNA.revcomp(seq));



    //testing
    $scope.x = {};
    $scope.x.obj = {};
    $scope.x.obj.y = "hey";
    //works
    console.log($scope.$eval('x.obj.y'));
    //works
    console.log($parse('x.obj.y')($scope));


    $scope.sequence = 'aaatttgggcccatgcta';

    //todo - if possible, assign these directly from DOM using a directive...
    $scope.$watch('sequence', function(newval, oldval) {
        $scope.rna = DNA.transcribe(newval);
    });

    $scope.$watch('rna', function (newval, oldval) {
        $scope.protein = DNA.translate(newval);
    });

    //Digest
    $scope.digestSeq = 'acaacgtctcacggatccagtcggaattctacatgcatcgatcgacggatccagatcgactagc';
    $scope.digestEnz = Digest.enzymes.EcoRI;
    $scope.digestFrags = Digest.digest($scope.digestSeq, $scope.digestEnz, false);

    $scope.$watch('digestEnz', function() {
        $scope.digestFrags = Digest.digest($scope.digestSeq, $scope.digestEnz, false);
    });


    //PCR predict

    $scope.pcr_demoSets = [
        {
            primers : ['tatcgatcgta', 'gatcgatcgat'],
            backbone : 'cccccccccccagctacgatcgataaaaaaaaaaattttttttttttgatcgatcgatagctaggggggggggggg'
        },
        {
            primers : ['CGGATCGATCGATCGT', 'GATCGATCGATAC'],
            backbone : 'tttttttttttttttttACGATCGATCGATCCGggggggggggggggggggcccccccccccccccGATCGATCGATACaaaaaaaaaaaaaa'
        },
        {
            primers : ['GTAGTAGTAGTAGTA', 'TCGATCGATCGATCGAA'],
            backbone : 'ttttttttttttttttttttttttttttttttttttttTACTACTACTACTACgggggggggggggggggggggggggggggggggggggggggggggggggggggggggccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccTCGATCGATCGATCGAaaaaaaaaaaaaaaaaaaaaa'
        },
        {
            primers : ['GTAGTAGTAA', 'TCGATCGAAA'],
            backbone : 'ttttttttttttttttttttttttttttttttttttttTACTACTACTACTACgggggggggggggggggggggggggggggggccccccccccccccccccccccccccccccTCGATCGATCGATCGAaaaaaaaaaaaaaaaaaaaaa'
        }
    ];

    $scope.setPCR = function(setInd) {
        var pcrSet = $scope.pcr_demoSets[setInd];
        $scope.primers = pcrSet.primers;
        $scope.backbone = pcrSet.backbone;
    };

    $scope.setPCR(2);
    
    
    
    console.log(PCR.findAnnealSimple($scope.backbone, $scope.primers[0]));


    //ligation

    //first 4 should all be the same
    $scope.ligate_demoSets = [
        ['aaaaaaaaaaA^CATG_', '^CATG_Tttggttggttgg'],
        ['aaaaaaaaaaA^CATG_', 'ccaaccaaccaaA^CATG_'],
        ['^CATG_Ttttttttttt', '^CATG_Tttggttggttgg'],
        ['^CATG_Ttttttttttt', 'ccaaccaaccaaA^CATG_'],
        ['aaaaaaaaaaA^CATG_', 'ggggggA^CATG_'],
        ['aaaaaaaaaaA^CATG_', 'gtcatcgatcagt_GTAC^'],
        ['aaaaaaaaaaA^CATG_T', 'A^CATG_Tacgatagcattaagcgt'],
        ['_gggg^tttttttttttt', 'aaaaa^cccc_'],
        ['aaaaaaaaaaaaaa|', '|ttggttggttgg']
    ];

    $scope.setLigate = function(setInd) {
        $scope.fragments = $scope.ligate_demoSets[setInd];
    };

    $scope.setLigate(0);


    $scope.$watch('fragments', function () {
        $scope.ligated = PCR.ligate($scope.fragments);
    }, true);


    $scope.clothoFunctions = [];
    Clotho.query({schema: "Function"}).then(function(result) {
        console.log(result);
        $scope.clothoFunctions = result;
    });


    //Construction Files

    $scope.constructionFileDemo = $http.get('/models/construction_demo.json').then(function(data) { return data.data });

    $scope.$watch('constructionFileDemo', function (newval) {
        if (!newval) return;
        $scope.constructionProductDemo = Construction.process(newval)
    });

    $scope.getConstructionFile = function(selection) {

        var files = [
            '/models/construction_demo.json',
            '/models/construction_gfp.json',
            '/models/construction_vio.json'
        ];

        $http.get(files[selection]).then(function(data) {
            $scope.constructionFile = data.data
        });
    };

    $scope.getConstructionFile(1);

    $scope.reprocess = function(file) {
        console.log('reprocess callback from controller!');
    };

    $scope.hideOpts = {};

}]);