'use strict';

$clotho.extensions.controller('clothoIntro_WritingCtrl', function($scope, $focus, $timeout, $http, DNA, Digest, PCR, Construction, Clotho) {

    $scope.getConstructionFile = function(selection) {

        var files = [
            'models/construction-old/construction_demo.json',
            'models/construction-old/construction_gfp.json',
            'models/construction-old/construction_vio.json',
            'models/construction-old/construction_skeleton.json'
        ];

        $http.get(files[selection]).then(function(data) {
            $scope.constructionFile = data.data
        });
    };

    $scope.reactionIndex = 1;

    $scope.$watch('reactionIndex', function() {
        $scope.getConstructionFile($scope.reactionIndex);
    });

    $scope.reactions = Construction.reactions;

});