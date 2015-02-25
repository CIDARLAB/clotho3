'use strict';

$clotho.extensions.controller('clothoIntro_SimulationCtrl', function($scope, $http, DNA, Digest, PCR, Construction, Clotho) {

    $http.get('models/construction-old/construction_gfp.json').then(function(data) {
        $scope.constructionFile = data.data
    });

    $scope.reactions = Construction.reactions;

});