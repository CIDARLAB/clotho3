'use strict';

describe('Controller: TestPlaylistimportCtrl', function () {

  // load the controller's module
  beforeEach(module('clothoApp'));

  var TestPlaylistimportCtrl,
    scope;

  // Initialize the controller and a mock scope
  beforeEach(inject(function ($controller, $rootScope) {
    scope = $rootScope.$new();
    TestPlaylistimportCtrl = $controller('TestPlaylistimportCtrl', {
      $scope: scope
    });
  }));

  it('should attach a list of awesomeThings to the scope', function () {
    expect(scope.awesomeThings.length).toBe(3);
  });
});
