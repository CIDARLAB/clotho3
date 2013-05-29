'use strict';

var testDynCtrl = /* Application.Controllers.controller('testDynCtrl', ['$scope', 'Clotho',  */
function($scope, Clotho) {

    $scope.editModel = function(uuid) {
        Clotho.edit(uuid);
    };

   $scope.hello = "world";
};

testDynCtrl.template = function(Clotho) {
 //return Clotho.get_url('show_template.html');
 return 'dynamic/dynamic-partial.html';
 };

dynamicCtrl.resolve = function(Clotho, $q, $timeout) {
    var resolved = {};
    var deferred = $q.defer();

    $q.all([
        resolved.model = Clotho.get('inst_first'),
        resolved.template_url = Clotho.get_url('show_template.html'),

        //testing - for exaggeration
        $timeout(function () {console.log('timeout')}, 1500)

    ]).then(function(values) {
        deferred.resolve(resolved);
    });

    return deferred.promise;
};