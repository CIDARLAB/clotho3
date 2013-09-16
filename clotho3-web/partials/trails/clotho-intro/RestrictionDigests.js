'use strict';

Application.Extensions.controller('clothoIntro_RestrictionDigestsCtrl', ['$scope', '$focus', '$timeout', '$dialog', 'Digest', 'Clotho', '$q', function($scope, $focus, $timeout, $dialog, Digest, Clotho, $q) {
    $scope.Digest = Digest;

    $scope.demo = {};
    $scope.demo.digestEnz = Digest.enzymes.EcoRI;
    $scope.demo.digestSeq = 'acaacgtctcacggatccagtcggaattctacatgcatcgatcgacggatccagatcgactagc';


    $scope.showMeHow = function () {
        $dialog.messageBox('Variables', 'The Search Bar creates a Rhino JavaScript environment in which you can run declare variables. We\'ll go ahead and declare <code>MySeq</code> and <code>MyEnz</code> for you behind the scenes.', [{label: "OK", cssClass: "btn-primary", result: true}]).open()
        .then(function() {
            var backgroundVars = $q.all([
                Clotho.submit('var mySeq = "acgatcgaatatACAGTGacgtactgatcga"'),
                //todo - update
                Clotho.submit('var myEnz = {match: "ACAGTG"}')
            ]);
            return $focus.typeOutSearch('clotho.run("digest", [mySeq, myEnz])')
        })
        .then();
    }
}]);