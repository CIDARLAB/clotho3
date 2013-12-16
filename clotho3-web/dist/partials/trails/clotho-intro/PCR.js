'use strict';

$clotho.extensions.controller('clothoIntro_PCRCtrl', ['$scope', '$focus', '$timeout', '$dialog', 'DNA', 'Digest', 'PCR', 'Clotho', 'Searchbar', function($scope, $focus, $timeout, $dialog, DNA, Digest, PCR, Clotho, Searchbar) {

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

        $dialog.messageBox('Your Turn!', 'We\'re leaving the variable declaration up to you this time. Try just copying some sequences from the PCR tool.', [{label: "No! I'm lazy.", cssClass: "", result: false}, {label: "OK", cssClass: "btn-primary", result: true}]).open()
        .then(function(result) {
            if (result)
                $('#clothoCommandBarInput').focus();
            else {
                Searchbar.submit('var myBackbone = "cccccccccccagctacgatcgataaaaaaaaaaatttttttttttgatcgatcgatagctaggggggggggggg"');
                Searchbar.submit("var myPrimers = ['tatcgatcgta', 'gatcgatcgat']");

                $dialog.messageBox('ಠ_ಠ', '<p>We\'re creating <code>myBackbone</code> and <code>myPrimers</code> as we speak!</p>', [{label: "Efficiency is intelligent laziness!", cssClass: "btn-primary", result: true}]).open()
                .then(function() {
                    var str = 'PCR.predict(myBackbone, myPrimers)';
                    $focus.typeOutSearch(str);
                    $('#clothoCommandBarInput').focus();
                })


            }
        });
    };
}]);