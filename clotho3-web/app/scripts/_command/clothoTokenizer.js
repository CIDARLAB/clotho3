angular.module('clotho.tokenizer')
/**
 * @ngdoc directive
 * @name clotho-tokenizer
 * @description
 * Wrapper for the clotho tokenizer. Creates token collection and autocomplete with dropdown
 * @example
 * <clotho-tokenizer ng-model="myModel" placeholder="Construct a command"></clotho-tokenizer>
 */
	.directive('clothoTokenizer', function ($parse, clothoTokenCollectionFactory, Debug) {

		var Debugger = new Debug('clothoTokenizer', '#ee7711');

		return {
			restrict: 'A',
			replace: true,
			require: 'ngModel', //avoid isolate scope so model propagates correctly
			templateUrl: "views/_command/tokenizer.html",
			link: function clothoTokenizerLink(scope, element, attrs, ngModelCtrl) {

				scope.placeholder = attrs.placeholder;

				var startingTags = $parse(attrs.startingTags)(scope);

				var completeQuery = '';

				scope.tokenCollection = new clothoTokenCollectionFactory(startingTags);

				/* updates + watches */

				function updateModel () {
					Debugger.log('updating model (query, tokens)', completeQuery, scope.tokenCollection.tokens);
					ngModelCtrl.$setViewValue({
						query: completeQuery,
						tokens : scope.tokenCollection.tokens
					});
				}

				function resetModel () {
					completeQuery = '';
					scope.tokenCollection.removeAll();
				}

				ngModelCtrl.$render = function () {
					//fixme - model is never actually updated, so we just reset it here assuming submit triggered $render
					resetModel();
					updateModel();
				};

				scope.$watchCollection('tokenCollection.tokens', function () {
					// if tokens changed, query should reflect it.
					// token is either just text, or a full sharable
					// can assume input element is irrelevant
					completeQuery = '';
					angular.forEach(scope.tokenCollection.tokens, function(token) {
						completeQuery  += token.readable() + ' ';
					});

					//update parent model
					updateModel();
				});

				/* functionality */

				scope.addToken = function (item) {
					Debugger.log('TOKENIZER_LINK adding token', item);
					scope.tokenCollection.addToken(item);
				};

				scope.removeToken = function (index, model) {
					Debugger.log('TOKENIZER_LINK removing token', index);
					scope.tokenCollection.removeToken(index);
				};

				scope.tokenActive = function (index) {
					return scope.tokenCollection.isActive(index);
				};

				scope.toggleTokenActive = function (index, event) {
					event.preventDefault();
					scope.tokenCollection.toggleActive(index);
				};

				scope.focusInput = function () {
					element[0].querySelector('[clotho-autocomplete]').focus();
				};
			}
		}
	});