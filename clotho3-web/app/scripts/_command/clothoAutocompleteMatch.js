angular.module('clotho.tokenizer')
	/*
	 * internal directive which displays an autocompletion
	 */
	.directive('clothoAutocompleteMatch', function (ClothoSchemas) {
		return {
			restrict:'EA',
			replace: true,
			scope:{
				index:'=',
				match:'=',
				//active : '@',
				query:'=',
				passedPlacement : '@'
			},
			templateUrl : 'views/_command/autocompleteMatch.html',
			link : function (scope, element, attrs) {
				scope.$watch(function () {
					return attrs['active'];
				}, function (newval) {
					scope.active = scope.$eval(newval);
				});

				scope.$watch('match', function (match) {
					scope.iconClass = ClothoSchemas.determineSharableIcon(ClothoSchemas.dirtyDetermineType(match));
				});
			}
		};
	})

	/*
	for autocomplete list, bold text matching query. requires that text bound is HTML not a string
	*/
	.filter('clothoAutocompleteHighlight', function() {

		function escapeRegexp(queryToEscape) {
			return queryToEscape.replace(/([.?*+^$[\]\\(){}|-])/g, "\\$1");
		}

		return function(matchItem, query) {
			return query ? matchItem.replace(new RegExp(escapeRegexp(query), 'gi'), '<strong>$&</strong>') : matchItem;
		};
	});