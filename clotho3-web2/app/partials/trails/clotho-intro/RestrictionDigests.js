'use strict';

angular.module('clotho.extensions').controller('clothoIntro_RestrictionDigestsCtrl', ['$scope', '$focus', '$timeout', '$dialog', 'Digest', 'Clotho', '$q', 'Searchbar', function($scope, $focus, $timeout, $dialog, Digest, Clotho, $q, Searchbar) {
    $scope.Digest = Digest;

    $scope.demo = {};
    $scope.demo.digestEnz = Digest.enzymes.EcoRI;
    $scope.demo.digestSeq = 'acaacgtctcacggatccagtcggaattctacatgcatcgatcgacggatccagatcgactagc';


    $dialog.messageBox('On to the Next Chapter...', '<p>Moving forward, we\'ll focus on using Clotho tools in sequence manipulation and other synthetic biology applications.</p> <p><b>We\'ll show you the tools, and click \"Show Me How\" to see how to use them in the Command Bar.</b></p>', [{label: "OK", cssClass: "btn-primary", result: true}]).open();


    $scope.showMeHow = function () {
        $dialog.messageBox('Variables', '<p>You can create variables in the Command Bar that Clotho will remember.</p><p>For example: <code>var mySeq = "acgatcgaatatGAATTCacgtactgatcga"</code></p><p>And using a Clotho-defined enzyme: <code>var myEnz = Digest.enzymes.EcoRI</code></p><p>We\'ll go ahead and declare <code>MySeq</code> and <code>MyEnz</code> for you behind the scenes.</p>', [{label: "OK", cssClass: "btn-primary", result: true}]).open()
        .then(function() {
            return $q.all([
                Searchbar.submit('var mySeq = "acgatcgaatatGAATTCacgtactgatcga"'),
                Searchbar.submit('var myEnz = Digest.enzymes.EcoRI'),
                $focus.typeOutSearch('Digest.digest(mySeq, myEnz)')
            ]);
        })
        .then(function() {
            $focus.focusSearch();
        });
    }
}]);