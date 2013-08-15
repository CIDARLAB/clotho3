'use strict';

/*This is meant for read-only modifications to the model - Can't use ngModel $formatters with promises (angular-1.1.5). Calling $setViewValue affects model and propagates, which is often undesired (e.g. if revcomp a sequence for display, don't want to change model).

NB if DO want to update model: ngModel weirdness in isolate scope -- call $parent.model
*/

//todo - broaden for other clotho function use - decide how to pass in
Application.Foundation.directive('clothoRun', ['Clotho', function(Clotho) {

    var inputsVal = {input: true, textarea : true, select: true};

    return {
        restrict : 'A',
        require : 'ngModel',
        scope : true,
        link: function (scope, element, attrs, ngModel) {

            //config
            var useVal = !!inputsVal[angular.lowercase(element[0].nodeName)];
            if (useVal) {
                //avoid flicker
                ngModel.$render = angular.noop;
            }

            //command, args
            scope.$watch(function() {
                return attrs.clothoRun
            }, function(newval, oldval) {
                console.log(newval);

                //todo - better handling

            });

            //model changes
            scope.$watch(function() {
                return ngModel.$modelValue
            }, function(newval, oldval) {
                console.log(newval);
                runFunction(newval);
            });


            var runFunction = function(input) {
                return Clotho.run(attrs.clothoRun, input).then(function(result) {
                    updateElement(result);
                });
            };

            var updateElement = function(newval) {
                var method = useVal ? 'val' : 'text';
                element[method](newval);
            };
        }
    }
}]);