'use strict';

$clotho.extensions.controller('clothoIntro_PCRCtrl', function($scope, $focus, $timeout, $modal, DNA, Digest, PCR, Clotho, CommandBar) {

    $scope.pcr_demoSets = [
        {
            primers : ['tatcgatcgta', 'gatcgatcgat'],
            backbone : 'ctcgatagcTACGATCGATAacgtagatcgatcagctgctagctagctactgaGATCGATCGATacgtacgtagctacgtacgac'
        },
        {
            primers : ['CGGATCGATCGATCGT', 'GATCGATCGATAC'],
            backbone : 'actgactacgatcgatACGATCGATCGATCCGatcgacgtaggacGATCGATCGATACacgatcgatcagctacgatcgatcgac'
        },
        {
            primers : ['GTAGTAGTAGTAGTA', 'TTCGATCGATCGATCGA'],
            backbone : 'acgtacgtacgTACTACTACTACTACagctagctacgatcgatcgatcgagcatctagctatcgtagctagcgactacgtacgtacgtcgatcgatcgatcgatcgatcgatcgactgatcgactgactgactagctacgatctacgtagctTCGATCGATCGATCGAtagctagctagctagctcga'
        },
        {
            primers : ['GTAGTAGTAA', 'TCGATCGAA'],
            backbone : 'gatcgatcgactgactagcatcgatcgatcgatcgatcgatcgagctacactgactactactattacTACTACTACTACTACacgtacgactacgacTCGATCGATCGATCGAactactacgatcgatcgatcgatcatc'
        }
    ];

    $scope.setPCR = function(setInd) {
        var pcrSet = $scope.pcr_demoSets[setInd];
        $scope.primers = pcrSet.primers;
        $scope.backbone = pcrSet.backbone;
    };

    $scope.setPCR(2);

    $scope.showMeHow = function() {

        $modal.messageBox('Your Turn!', 'We\'re leaving the variable declaration up to you this time. Try just copying some sequences from the PCR tool.', [{label: "No! I'm lazy.", cssClass: "", result: false}, {label: "OK", cssClass: "btn-primary", result: true}])
        .result
        .then(function(result) {
            if (result)
                $('#clothoCommandBarInput').focus();
            else {
                CommandBar.submit('var myBackbone = "cccccccccccagctacgatcgataaaaaaaaaaatttttttttttgatcgatcgatagctaggggggggggggg"');
	            CommandBar.submit("var myPrimers = ['tatcgatcgta', 'gatcgatcgat']");

                $modal.messageBox('ಠ_ಠ', '<p>We\'re creating <code>myBackbone</code> and <code>myPrimers</code> as we speak!</p>', [{label: "Efficiency is intelligent laziness!", cssClass: "btn-primary", result: true}])
                .result
                .then(function() {
                    var str = 'PCR.predict(myBackbone, myPrimers)';
                    $focus.typeOutSearch(str);
                    $('#clothoCommandBarInput').focus();
                })


            }
        });
    };
});