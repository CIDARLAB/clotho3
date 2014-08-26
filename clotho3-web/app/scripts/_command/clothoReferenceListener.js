/**
 * @ngdoc directive
 * @name clothoReferenceListener
 * @description
 * Listens for keydown of delimiter symbol (@), and triggers a callback when typed
 *
 * @attr clothoReferencePreventDefault
 * @attr clothoReferenceCallback
 * Passed $event, $element, $delimiter
 */
angular.module('clotho.tokenizer')
	.directive('clothoReferenceListener', function (ClothoReferenceDelimiter, $parse) {
		return function (scope, element, attrs) {
			element.on('keydown', function (evt) {
				if (event.which === ClothoReferenceDelimiter.keycode) {
					if (scope.$eval(attrs.clothoReferencePreventDefault)){
						evt.preventDefault();
					}

					$parse(attrs.clothoReferenceCallback)(scope, {
						$event : evt,
						$element : element,
						$delimiter : ClothoReferenceDelimiter.symbol
					});
				}
			});
		};
	});