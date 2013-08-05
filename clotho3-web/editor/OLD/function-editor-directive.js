'use strict';

Application.Editor.directive('functionEditor', ['Clotho', '$compile', '$parse', function(Clotho, $compile, $parse) {

    return {
        restrict: 'CA',
        require: '^form',
        templateUrl: '/editor/function-partial.html',
        scope: {
            uuid: '='
        },
        controller: function($scope, $element, $attrs) {

            //if we're not linking, just pull it from the attrs
            if (typeof $scope.uuid == 'undefined') {
                console.log("id not in controller, pulling");
                $scope.uuid = $attrs.uuid;
            }

            $scope.editing = $scope.editing || {};

            $scope.editMode = false;
            $scope.formDirty = false;

            $scope.languages = ['javascript', 'python', 'java'];

            $scope.paramTypes = [
                {name:'Object', type:'Type'},
                {name:'String', type:'Type'},
                {name:'Integer', type:'Type'},
                {name:'Sequence', type:'Schema'},
                {name:'Person', type:'Schema'},
                {name:'Institution', type:'Schema'}
            ];




            function emptyParam() {return {"type" : "", "name" : "", "test" : {"uuid" : ""}}}

            $scope.addParam = function() {
                if (angular.isEmpty($scope.editing.params)) {$scope.editing.params = [];}
                $scope.editing.params.push(emptyParam());
            };



        },
        link: function link(scope, element, attrs, ngForm) {

            //e.g. scope.formConst.$setPristine()
            scope.formConst = $parse(attrs.name)(scope);

            //switch to 'edit' mode
            scope.edit = function() {
                scope.editMode = true;
            };

            //discard edits
            scope.reset = function() {
                scope.formConst.$setPristine();
                Clotho.get(scope.uuid).then(function(result) {
                    scope.function = result;
                });
            };

            //save edits, switch to 'view'
            scope.save = function() {
                Clotho.set(scope.function);
                scope.editMode = false;
            };

            //discard edits, switch to 'view'
            scope.discard = function() {
                scope.reset();
                scope.editMode = false;
            };
        }
    }
}]);