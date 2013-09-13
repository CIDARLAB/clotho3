'use strict';

Application.Extensions.controller('writingTrailsCtrl', ['$compile', '$scope', 'Clotho', '$http', function($compile, $scope, Clotho, $http) {
    $scope.getSchema = function(schemaName) {
        return Clotho.query({schema: "ClothoSchema", name: schemaName})
            .then(function(results) {
                console.log(results[0]);
                return results[0]
            });
    };

    $scope.getTrail = function(trailName) {
        return Clotho.query({schema: 'Trail', name: trailName})
            .then(function(results) {
                console.log(results[0]);
                return results[0]
            });
    };

    //temporary hack while Trail schema can't be imported
    $scope.getTrailHack = function() {
        return $http.get('/models/_NOT TO SYNC/definition_trail.json').then(function(data) {
            return data.data;
        });
    }


}]);

Application.Extensions.controller('writingTrailsUICtrl', ['Clotho', '$scope', '$focus', '$dialog', function(Clotho, $scope, $focus, $dialog) {
    $scope.videoDemo = {
            id : '<11 char videoID>',
            params : '<youtubeParamsObj or URL>',
            autoplay : 'true|false',
            mini : 'true|false'
        };

        $scope.showSearchTypeOut = function () {
            $focus.typeOutSearch('Demo search text');
        };

        $scope.showDialog = function() {
            $dialog.messageBox('Demo Dialog', 'This is some simple text. It can contain <code>HTML</code>', [{label: "OK", cssClass: "btn-primary", result: true}]).open()
                .then(function() {
                    return $dialog.messageBox('Chaining', 'The $dialog and $focus services return promises so you can chain together UI elements based on their completion / response', [{label: "OK", cssClass: "btn-primary", result: true}]).open()
                })
                .then(function() {
                    return $focus.typeOutSearch("$focus demo")
                })

        };

        //demo video modal
        $scope.openVideo = function() {
            $dialog.video('ivif214cQLE', {}).open();
        };

        $scope.schemas = [];
        Clotho.query({"schema": "Schema"}).then(function(data) {
            $scope.schemas = data;
        });


        $scope.dropdownItems = [
            "The first choice!",
            "And another choice for you.",
            "but wait! A third!"
        ];


}]);


Application.Extensions.controller('writingTrailsClothoCtrl', ['Clotho', '$scope', '$focus', '$dialog', 'Digest', function(Clotho, $scope, $focus, $dialog, Digest) {

    $scope.Digest = Digest;
    //Digest
    $scope.digestSeq = 'acaacgtctcacggatccagtcggaattctacatgcatcgatcgacggatccagatcgactagc';
    $scope.digestEnz = Digest.enzymes.EcoRI;


}]);