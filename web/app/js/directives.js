'use strict';

/* Directives */

angular.module("editorDirectives", [])
    .directive('form', function () {
        return {
            restrict: 'E',
            require: 'form',
            link: function (scope, element, attrs, ctrl) {
                //not really used yet but here for example
                scope.$watch(attrs.name + '.$dirty', function (newValue, oldValue) {
                    if (newValue != oldValue && newValue === true) {
                        scope.$emit('formDirty', attrs.name);
                    }
                });
            }
        }
    });