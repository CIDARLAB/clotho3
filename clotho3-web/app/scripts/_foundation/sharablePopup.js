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
 * todo - need to persist type rather than run() the check every time popup opens
 *
 */
	.directive('sharablePopup', function ($clothoPopup) {
		return $clothoPopup('sharablePopup');
	})
	.directive('sharablePopupInner', function (Clotho, ClothoSchemas, $timeout, $injector) {

		var editorPresent = $injector.has('clothoEditorDirective');

		//future - check for $route and /executor --- if not present then executor should popup a modal

		return {
			restrict: 'EA',
			replace: true,
			scope: {
				sharableId: '=?',
				sharableModel : '=?',
				placement: '@popupPlacement',
				reposition : '&'
			},
			templateUrl: 'views/_foundation/sharableBasicFieldsPopup.html',
			link : function (scope, element, attrs, nullCtrl) {

				//given a model (should be pruned), sets related scope variables
				function setSharable (model) {
					scope.sharable = model;

					//let's make sure this is right and defer to the server version
					ClothoSchemas.determineSharableType(model)
					.then(function (detType) {
						scope.type = detType;
					}, function (err) {
						scope.type = ClothoSchemas.dirtyDetermineType(model);
					})
					.then(function () {
						scope.iconClass = ClothoSchemas.determineSharableIcon(scope.type);
						scope.labelClass = 'label-' + ClothoSchemas.typeToColorClass(scope.type);
					});


					if (ClothoSchemas.isSchema(model)) {
						scope.isSchema = true;
						//must have a proper schema to download schema dependencies
						Clotho.get(model.id, {mute: true})
						.then(function (fullModel) {
							ClothoSchemas.getSuperclassFields(fullModel)
							.then(function (otherFields) {
								scope.schema = fullModel;
								scope.inheritedFields = otherFields;
							});
						});
					}
				}

				scope.editorPresent = editorPresent;

				scope.$watch('sharableModel', function ( val, oldval ) {
					if (!!val) {
						setSharable(val);
						//no need to reposition if we pass in the model, since the initial $digest will fill it

						if (val.id) {
							Clotho.get(val.id, {mute : true}).then(function (retrievedSharable) {
								scope.fullSharable = retrievedSharable;
							});
						}
					}
				});

				scope.$watch('sharableId', function ( val, oldval ) {
					if (!!val) {
						Clotho.get(val, {mute : true}).then(function (retrievedSharable) {
							if (!angular.isEmpty(retrievedSharable)) {
								scope.fullSharable = retrievedSharable;
								setSharable( ClothoSchemas.pruneToBasicFields(retrievedSharable) );
								//if we're getting it remotely, the size probably changed, so let's reposition it.
								scope.reposition();
							}
						});
					}
				});

				scope.showView = function (evt) {
					scope.activeView = scope.activeView ? '' : scope.type;
					$timeout(function () {
						scope.reposition();
					});
				};

				scope.edit = Clotho.edit;

				scope.$on('$destroy', function () {
          //todo - persist state (type) in localstorage or something
				});
			}
		};
	});
