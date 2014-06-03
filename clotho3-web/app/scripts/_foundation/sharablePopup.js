angular.module('clotho.clothoDirectives')
/**
 * @ngdoc directive
 * @name sharable-popup
 *
 * @description Displays a popup showing the basic fields of an instance. Appended to Body. Either pass Model or ID of sharable, model gets priority (should only use pruned fields or may be very large)
 *
 * @note - need to create isolate scope so properties don't overlap
 *
 * @attrs
 * sharablePopupModel
 * sharablePopupId
 * sharablePopupPlacement
 * sharablePopupTrigger (none | click | mouseenter | focus)
 * sharablePopupOpen
 *
 *
 */
	.directive('sharablePopup', function ($clothoPopup) {
		return $clothoPopup('sharablePopup');
	})
	.directive('sharablePopupInner', function (Clotho, ClothoSchemas, $injector) {

		var editorPresent = $injector.has('clothoEditorDirective');

		return {
			restrict: 'EA',
			replace: true,
			scope: {
				sharableId: '=?',
				sharableModel : '=?',
				placement: '@',
				reposition : '&'
			},
			templateUrl: 'views/_foundation/sharableBasicFieldsPopup.html',
			link : function (scope, element, attrs, nullCtrl) {

				//given a model (should be pruned), sets related scope variables
				function setSharable (model) {
					scope.sharable = model;
					scope.type = ClothoSchemas.determineSharableType(model);
					scope.iconClass = ClothoSchemas.determineSharableIcon(scope.type);
					scope.labelClass = 'label-' + ClothoSchemas.typeToColorClass(scope.type);

					if (ClothoSchemas.isSchema(model)) {
						scope.isSchema = true;
						//must have a proper schema to download schema dependencies
						Clotho.get(model.id, {mute: true})
							.then(function (fullModel) {
								ClothoSchemas.downloadSchemaDependencies(fullModel)
									.then(function (finalSchema) {
										scope.schema = finalSchema;
									});
							})

					}
				}

				scope.editorPresent = editorPresent;

				scope.$watch('sharableModel', function ( val, oldval ) {
					if (!!val) {
						setSharable(val);
						//no need to reposition if we pass in the model, since the initial $digest will fill it
					}
				});

				scope.$watch('sharableId', function ( val, oldval ) {
					if (!!val && angular.isEmpty(scope.sharableModel)) {
						Clotho.get(val, {mute : true}).then(function (retrievedSharable) {
							scope.fullSharable = retrievedSharable;
							setSharable( ClothoSchemas.pruneToBasicFields(retrievedSharable) );
							//if we're getting it remotely, the size probably changed, so let's reposition it.
							scope.reposition();
						});
					}
				});

				//remove focus from autocomplete input... hard to hide popup.
				scope.toggleSchema = function (evt) {
					evt.preventDefault();
					scope.showingSchema = !scope.showingSchema;
					scope.reposition();
				};

				scope.edit = Clotho.edit;

				scope.$on('$destroy', function () {
				})
			}
		};
	});